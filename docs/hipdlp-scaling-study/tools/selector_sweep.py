#!/usr/bin/env python3
"""Fair CPU iteration A/B, with immutable sources copied into scratch.

Names must be supplied explicitly. Held-out mode rejects study-list names and
known filesystem aliases. Each run records the command, complete log, and JSON.
Time-limited pairs are censored in the report, never treated as iteration wins.
"""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adapter', type=Path, required=True)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--source', type=Path, default=Path('/srv/mps_da'))
    parser.add_argument('--names', type=Path, required=True)
    parser.add_argument('--held-out', action='store_true')
    parser.add_argument('--iterations', type=int, default=500000)
    parser.add_argument('--time-limit', type=float, default=600)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--configs', default='A,L,S')
    args = parser.parse_args()
    names = args.names.read_text().split()
    aliases = {'neos-5052403': 'neos5052403',
               'neos-5251015': 'neos5251015',
               'physiciansched3-3': 'psched3-3'}
    study = Path(__file__).resolve().parents[1]
    training = set((study / 'data/instances.txt').read_text().split())
    if args.held_out:
        overlap = set(names) & (training | set(aliases.values()))
        if overlap:
            parser.error(f'held-out names overlap training: {sorted(overlap)}')
    lp_dir = args.scratch / 'lp'
    logs = args.scratch / 'runs'
    lp_dir.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name in names:
        source = args.source / (aliases.get(name, name) + '.mps.gz')
        copied = lp_dir / (name + '.mps.gz')
        shutil.copyfile(source, copied)
        with copied.open('rb') as stream:
            compressed = stream.read(2) == b'\x1f\x8b'
        with (gzip.open(copied, 'rb') if compressed else copied.open('rb')) as src:
            with (lp_dir / (name + '.mps')).open('wb') as dest:
                shutil.copyfileobj(src, dest)
        manifest.append({'name': name, 'source': str(source),
                         'sha256': hashlib.file_digest(copied.open('rb'), 'sha256').hexdigest()})
    (args.scratch / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

    def run(item):
        name, config = item
        prefix = logs / f'{name}__{config}'
        command = [str(args.adapter.resolve()), '--instance', str(lp_dir / (name + '.mps')),
                   '--backend', 'cpu', '--scaling', 'legacy' if config == 'L' else 'auto',
                   '--step-size', 'fixed', '--presolve', 'off', '--kkt-tolerance', '1e-4',
                   '--iteration-limit', str(args.iterations), '--time-limit', str(args.time_limit), '--log']
        if config == 'S':
            command.append('--scaling-selector')
        Path(str(prefix) + '.command.json').write_text(json.dumps(command) + '\n')
        with Path(str(prefix) + '.log').open('w') as output:
            process = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT)
        raw = Path(str(prefix) + '.log').read_text()
        start = raw.find('{\n  "schema"')
        result = json.loads(raw[start:]) if start >= 0 else {'status': 'ADAPTER_ERROR'}
        result.update(name=name, config=config, returncode=process.returncode,
                      iteration_budget=args.iterations, held_out=args.held_out)
        Path(str(prefix) + '.json').write_text(json.dumps(result, indent=2) + '\n')
        print(name, config, result['status'], result.get('iterations'), flush=True)
        return result

    items = [(name, config) for name in names for config in args.configs.split(',')]
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(run, items))
    (args.scratch / 'results.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
