#!/usr/bin/env python3
"""Copy benchmarks to scratch; resume CPU probes or serial CUDA comparisons."""
import argparse
import concurrent.futures
import json
from pathlib import Path
import shutil
import subprocess
from prepare_large import decode, digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--driver', type=Path, required=True)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--names', type=Path, required=True)
    p.add_argument('--source', type=Path, default=Path('/srv/mps_da'))
    p.add_argument('--backend', choices=['cpu', 'cuda'], default='cpu')
    p.add_argument('--configs', default='F')
    p.add_argument('--seconds', type=float, default=60)
    p.add_argument('--iterations', type=int, default=200)
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--model', type=Path)
    a = p.parse_args()
    if a.backend == 'cuda' and a.jobs != 1:
        p.error('CUDA timings require --jobs 1')
    names = a.names.read_text().split()
    lp = a.scratch / 'lp'
    runs = a.scratch / 'runs'
    lp.mkdir(parents=True, exist_ok=True)
    runs.mkdir(parents=True, exist_ok=True)
    aliases = {'neos-5052403': 'neos5052403', 'neos-5251015': 'neos5251015',
               'physiciansched3-3': 'psched3-3'}
    manifest_path = a.scratch / 'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else []
    by_name = {r['name']: r for r in manifest}
    for name in names:
        src = a.source / (aliases.get(name, name) + '.mps.gz')
        target = lp / (name + '.mps')
        if name in by_name and target.exists():
            if digest(src) != by_name[name]['sha256'] or digest(target) != by_name[name]['mps_sha256']:
                raise ValueError(f'input changed: {name}')
            continue
        copy = lp / (name + '.mps.gz')
        shutil.copyfile(src, copy)
        members, trailing = decode(copy, target)
        row = dict(name=name, source=str(src), sha256=digest(copy),
                   mps_sha256=digest(target), gzip_members=members,
                   non_gzip_trailing_bytes=trailing)
        manifest.append(row)
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    model_hash = digest(a.model) if a.model else None
    driver_hash = digest(a.driver)

    def run(item):
        name, config = item
        out = runs / f'{name}__{config}.json'
        cmd = [str(a.driver.resolve()), str(lp / (name + '.mps')), name,
               a.backend, config, str(a.seconds), str(a.iterations),
               str(a.model.resolve()) if a.model else '-']
        if out.exists():
            row = json.loads(out.read_text())
            if row['command'] != cmd or row['driver_sha256'] != driver_hash or row['model_sha256'] != model_hash:
                raise ValueError(f'run metadata changed: {out}')
            return row
        log = out.with_suffix('.log')
        with log.open('w') as stream:
            r = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT,
                               timeout=a.seconds + 180)
        raw = log.read_text()
        lines = [s[len('ONE_RESULT '):] for s in raw.splitlines() if s.startswith('ONE_RESULT ')]
        if r.returncode or len(lines) != 1:
            raise RuntimeError(f'driver failed: {log}')
        row = json.loads(lines[0])
        row.update(command=cmd, driver_sha256=driver_hash, model_sha256=model_hash)
        out.write_text(json.dumps(row, indent=2, allow_nan=False) + '\n')
        print(name, config, row['status'], row['iterations'],
              round(row['solver_wall_sec'], 4), row['chosen'], flush=True)
        return row
    configs = a.configs.split(',')
    # Rotate order across instances; keep all GPU solves serial.
    items = [(name, config) for i, name in enumerate(names)
             for config in configs[i % len(configs):] + configs[:i % len(configs)]]
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
        rows = list(pool.map(run, items))
    (a.scratch / 'results.json').write_text(json.dumps(rows, indent=2) + '\n')
    for row in manifest:
        if digest(Path(row['source'])) != row['sha256']:
            raise ValueError(f'source changed: {row["source"]}')
    (a.scratch / 'integrity.json').write_text(json.dumps(dict(
        source_hashes_unchanged=True, count=len(manifest), driver_sha256=driver_hash,
        model_sha256=model_hash), indent=2) + '\n')


if __name__ == '__main__':
    main()
