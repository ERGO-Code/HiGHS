#!/usr/bin/env python3
"""Filter matrix roundoff and declare bound-only columns without changing bounds/costs.

For these free-format NESO MPS files, keep all rows and all bound definitions.
Explicit zero objective declarations stop readers from silently omitting columns
that appear only in BOUNDS. Removing coefficients is an explicit model change.
"""
import argparse
import hashlib
import json
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def canonicalize(source, target, threshold):
    columns, bound_columns, objective = set(), [], None
    section = ''
    with source.open() as f:
        for line in f:
            tokens = line.split()
            if not tokens or line.startswith('*'):
                continue
            if not line[0].isspace():
                section = tokens[0]
                continue
            if section == 'ROWS' and tokens[0] == 'N' and objective is None:
                objective = tokens[1]
            elif section == 'COLUMNS':
                assert len(tokens) in (3, 5), tokens
                columns.add(tokens[0])
            elif section == 'BOUNDS':
                assert len(tokens) in (3, 4), tokens
                if tokens[2] not in columns:
                    columns.add(tokens[2])
                    bound_columns.append(tokens[2])
    assert objective is not None
    removed = 0
    section = ''
    emitted_columns = set()
    with source.open() as f, target.open('w') as out:
        for line in f:
            tokens = line.split()
            if not tokens or line.startswith('*'):
                out.write(line)
                continue
            if not line[0].isspace():
                if section == 'COLUMNS':
                    for name in bound_columns:
                        out.write(f'    {name}  {objective}  0\n')
                section = tokens[0]
                out.write(line)
                continue
            if section != 'COLUMNS':
                out.write(line)
                continue
            pairs = []
            for k in range(1, len(tokens), 2):
                if tokens[k] == objective or abs(float(tokens[k + 1])) > threshold:
                    pairs.extend(tokens[k:k + 2])
                else:
                    removed += 1
            if pairs:
                out.write('    ' + tokens[0] + '  ' + '  '.join(pairs) + '\n')
                emitted_columns.add(tokens[0])
            elif tokens[0] not in emitted_columns:
                # Keep the variable declaration even if every matrix entry was dropped.
                out.write(f'    {tokens[0]}  {objective}  0\n')
                emitted_columns.add(tokens[0])
    return dict(source=str(source), target=str(target), source_sha256=sha256(source),
                target_sha256=sha256(target), threshold=threshold, removed=removed,
                bound_only_columns=len(bound_columns), declared_columns=len(columns))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('source', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--threshold', type=float, default=1e-9)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for name in ('neso-2245', 'neso-2005'):
        record = canonicalize(args.source / (name + '.mps'), args.output / (name + '-canonical.mps'), args.threshold)
        print(json.dumps(record), flush=True)
        records.append(record)
        (args.output / 'manifest.json').write_text(json.dumps(records, indent=2) + '\n')


if __name__ == '__main__':
    main()
