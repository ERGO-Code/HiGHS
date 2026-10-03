#!/usr/bin/env python3
"""Analyze stair_diag's simplex-reference operators and saved pilot solutions.

Uses NumPy only. The active-subproblem spectrum is diagnostic evidence, not a
claim that a particular singular value determines the nonlinear restart rate.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def read_operator(path):
    with path.open() as stream:
        m, n, entries, weight = stream.readline().split()
        m, n, entries = int(m), int(n), int(entries)
        cols = np.array([[float(x) for x in stream.readline().split()] for _ in range(n)])
        rows = np.array([[float(x) for x in stream.readline().split()] for _ in range(m)])
        matrix = np.zeros((m, n))
        for _ in range(entries):
            row, col, value = stream.readline().split()
            matrix[int(row), int(col)] = float(value)
    return matrix, cols, rows, float(weight)


def spectrum(matrix):
    values = np.linalg.svd(matrix, compute_uv=False)
    nonzero = values[values > values[0] * 1e-10]
    return dict(shape=list(matrix.shape), rank=len(nonzero), sigma_max=float(values[0]),
                sigma_min_positive=float(nonzero[-1]), condition=float(nonzero[0] / nonzero[-1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-prefix', type=Path, required=True)
    parser.add_argument('--trajectories', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = {}
    for config in ('A', 'L'):
        matrix, cols, rows, weight = read_operator(Path(str(args.reference_prefix) + '-' + config + '.txt'))
        x, lower, upper, rc, col_scale, cost = cols.T
        y, rhs, row_upper, equality, row_scale, original_row, row_sign = rows.T
        active = (equality > 0) | (abs(y) > 1e-7)
        free = ((~np.isfinite(lower)) | (x > lower + 1e-7)) & ((~np.isfinite(upper)) | (x < upper - 1e-7))
        initial = np.clip(np.zeros_like(x), lower, upper)
        primal_term = weight * np.dot(x - initial, x - initial) / 2
        dual_term = np.dot(y, y) / (2 * weight)
        full_spectrum = spectrum(matrix)
        # Include weakly active columns with zero reduced cost. This avoids
        # assigning the reference basis's arbitrary degenerate bounds a unique
        # derivative, and provides a second active-set diagnostic alongside
        # the strictly interior-column spectrum.
        zero_rc_matrix = matrix[np.ix_(active, abs(rc) < 1e-8)]
        left_modes, _, _ = np.linalg.svd(zero_rc_matrix, full_matrices=False)
        values = dict(full_spectrum=full_spectrum,
                      active_free_spectrum=spectrum(matrix[np.ix_(active, free)]),
                      active_zero_rc_spectrum=spectrum(zero_rc_matrix),
                      weight=weight, primal_term=float(primal_term), dual_term=float(dual_term),
                      reference_Kw=float(full_spectrum['sigma_max'] * (primal_term + dual_term)))
        # stair has no synthetic slacks, so original columns equal transformed
        # columns. Check rather than silently using this on another model.
        trajectories = []
        for pilot_config in ('A', 'L'):
            for length in (200, 1000, 4000, 16000, 64000, 128000, 256000):
                path = args.trajectories / f'{pilot_config}-{length}.sol'
                if not path.exists():
                    continue
                with path.open() as stream:
                    nc, nr = map(int, stream.readline().split())
                    assert nc == len(x) and nr == len(y), 'only stair without slacks is supported'
                    data = np.array([float(line) for line in stream])
                xp = data[:nc] / col_scale
                yp = data[2 * nc:][original_row.astype(int)] * row_sign / row_scale
                pt = weight * np.dot(xp - initial, xp - initial) / 2
                dt = np.dot(yp, yp) / (2 * weight)
                error_energy = (weight * np.dot(xp - x, xp - x) + np.dot(yp - y, yp - y) / weight) / 2
                dual_error = (yp - y)[active]
                coefficients = left_modes.T @ dual_error
                dual_error_squared = float(np.dot(dual_error, dual_error))
                trajectories.append(dict(pilot=pilot_config, length=length, Dx=float(np.linalg.norm(xp - initial)),
                                         Dy=float(np.linalg.norm(yp)), primal_term=float(pt), dual_term=float(dt),
                                         Kw=float(full_spectrum['sigma_max'] * (pt + dt)),
                                         reference_error_energy=float(error_energy),
                                         slowest_dual_mode_fraction=float(coefficients[-1] ** 2 / dual_error_squared)
                                         if dual_error_squared else 0.0))
        values['trajectories'] = trajectories
        result[config] = values
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    for config, values in result.items():
        print(config, 'Kw', values['reference_Kw'], 'active spectrum', values['active_free_spectrum'])


if __name__ == '__main__':
    main()
