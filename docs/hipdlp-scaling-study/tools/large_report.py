#!/usr/bin/env python3
"""Summarize verified solves and fair elapsed times without rewarding failures."""
import argparse
import csv
import json
import math
from pathlib import Path


LABELS = {'A': 'auto (guard + CR/Ruiz/PC)', 'L': 'Ruiz+PC',
          'S': 'two-pilot selector', 'P': 'PC-only', 'R': 'Ruiz-only',
          'C': 'CR-only', 'O': 'off', 'AP': 'CR+PC', 'AR': 'CR+Ruiz'}


def verified(row):
    if row['status'] == 'OPTIMAL':
        return row['kkt_pass']
    if row['status'] in ('INFEASIBLE', 'UNBOUNDED', 'UNBOUNDED_OR_INFEASIBLE'):
        return row['certificate_pass']
    return False


def gm(values):
    return math.exp(sum(math.log(max(v, 1e-12)) for v in values) / len(values)) if values else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    rows = []
    for path in args.directory.rglob('*__*.json'):
        row = json.loads(path.read_text())
        if isinstance(row, dict) and 'config' in row:
            rows.append(row)
    configs = sorted({r['config'] for r in rows})
    names = sorted({r['name'] for r in rows})
    table = {(r['name'], r['config']): r for r in rows}
    if len(table) != len(rows) or len({r['backend'] for r in rows}) > 1:
        parser.error('choose one experiment directory; duplicate runs or mixed backends found')
    common = [n for n in names if all((n, c) in table and verified(table[n, c]) for c in configs)]
    read_total = sum(next(r['model']['read_wall_sec'] for r in rows
                         if r['name'] == name) for name in names)
    summary = {'measured_runs': len(rows), 'instances': len(names),
               'shared_input_read_total_sec': read_total,
               'common_verified_names': common, 'configurations': {}}
    lines = [f'{len(rows)} completed runs across {len(names)} instances.', '',
             '| Configuration | Verified | Measured total (s) | PAR2 (s) | Common solved total (s) |',
             '|---|---:|---:|---:|---:|']
    for config in configs:
        subset = [r for r in rows if r['config'] == config]
        passed = [r for r in subset if verified(r)]
        totals = {'completed': len(subset), 'verified': len(passed),
                  'verified_solved_names': sorted(r['name'] for r in passed),
                  'verified_solved_total_sec': sum(r['solver_wall_sec'] for r in passed),
                  'observed_total_sec': sum(r['solver_wall_sec'] for r in subset),
                  'par2_sec': sum(r['solver_wall_sec'] if verified(r) else 2 * r['time_limit'] for r in subset),
                  'common_solved_total_sec': sum(table[n, config]['solver_wall_sec'] for n in common),
                  'common_solved_time_geomean_sec': gm([table[n, config]['solver_wall_sec'] for n in common])}
        totals['observed_total_with_shared_input_read_sec'] = totals['observed_total_sec'] + read_total
        if subset and subset[0]['backend'] == 'cpu':
            totals['common_solved_iteration_geomean'] = gm([table[n, config]['iterations'] for n in common])
        against_a = []
        for r in subset:
            a = table.get((r['name'], 'A'))
            if a and verified(a) and verified(r):
                against_a.append({'name': r['name'], 'time_ratio_to_auto': r['solver_wall_sec'] / a['solver_wall_sec'],
                                  'seconds': r['solver_wall_sec'], 'auto_seconds': a['solver_wall_sec']})
        totals['time_ratio_geomean_to_auto'] = gm([r['time_ratio_to_auto'] for r in against_a])
        totals['worst_time_regressions_to_auto'] = sorted(against_a, key=lambda r: -r['time_ratio_to_auto'])[:5]
        totals['lost_solves_to_auto'] = [r['name'] for r in subset if not verified(r) and
            (r['name'], 'A') in table and verified(table[r['name'], 'A'])]
        totals['censored_time_regressions_to_auto'] = sorted([
            {'name': r['name'], 'status': r['status'],
             'time_limit_ratio_to_auto': r['time_limit'] / table[r['name'], 'A']['solver_wall_sec']}
            for r in subset if r['status'] == 'TIME_LIMIT' and
            (r['name'], 'A') in table and verified(table[r['name'], 'A'])],
            key=lambda r: -r['time_limit_ratio_to_auto'])
        totals['gained_solves_to_auto'] = [r['name'] for r in subset if verified(r) and
            (r['name'], 'A') in table and not verified(table[r['name'], 'A'])]
        summary['configurations'][config] = totals
        lines.append(f'| {LABELS[config]} | {len(passed)}/{len(subset)} | '
                     f'{totals["observed_total_sec"]:.3f} | {totals["par2_sec"]:.3f} | '
                     f'{totals["common_solved_total_sec"]:.3f} |')
    lines += ['', f'Common verified subset: {len(common)} instances: {", ".join(common)}.', '',
              'Measured total includes failed runs at their observed duration. PAR2 charges',
              'each failure twice its solver time limit. Input I/O and backend warmup are',
              'excluded; scaling, pilot work and the solver final KKT check are included.', '']
    (args.directory / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (args.directory / 'results.json').write_text(json.dumps(
        sorted(rows, key=lambda r: (r['name'], r['config'])), indent=2) + '\n')
    (args.directory / 'summary.md').write_text('\n'.join(lines))
    fields = ['name', 'config', 'backend', 'status', 'kkt_pass', 'certificate_pass',
              'solver_wall_sec', 'verification_wall_sec', 'iterations', 'guard_rejected',
              'time_limit', 'iteration_limit']
    with (args.directory / 'results.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r['name'], r['config'])))
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
