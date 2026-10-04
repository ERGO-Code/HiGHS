#!/usr/bin/env python3
"""Archive compact experiment evidence without LPs or large solution vectors."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('scratch', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name in ('verified-results.json', 'canonical-checked/manifest.json'):
        shutil.copyfile(args.scratch / name, args.output / Path(name).name)
    provenance = args.scratch / 'input-provenance.json'
    if provenance.exists():
        shutil.copyfile(provenance, args.output / provenance.name)
    repos = dict(hipdlp2=Path(__file__).resolve().parents[3],
                 hpr_julia=args.scratch / 'HPR-LP-current',
                 hpr_c=args.scratch / 'HPR-LP-C',
                 cupdlpx=args.scratch / 'cuPDLPx-current')
    metadata = dict(revisions={}, binaries={},
                    gpu='NVIDIA GeForce RTX 4070 Ti SUPER (16 GB)',
                    compiler_toolkit='CUDA 12.8',
                    acceptance='hipdlp-original-kkt/2 at 1e-4',
                    presolve='off', numerical_source_changes=False)
    for name, repo in repos.items():
        result = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'], capture_output=True, text=True, check=True)
        metadata['revisions'][name] = result.stdout.strip()
        if name in ('hpr_c', 'cupdlpx'):
            patch = subprocess.run(['git', '-C', str(repo), 'diff'], capture_output=True, text=True, check=True)
            (args.output / (name + '-export.patch')).write_text(patch.stdout)
    for relative in ('hpr-build/solve_mps_file', 'cupdlpx-current-build/cupdlpx', 'rewrite-adapter', 'neso-kkt', 'rewrite-profile'):
        path = args.scratch / relative
        metadata['binaries'][relative] = sha(path)
    (args.output / 'build-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    for file in (args.scratch / 'runs').glob('*/*/*'):
        if file.suffix == '.json' or file.name.endswith('_summary.txt'):
            destination = args.output / file.relative_to(args.scratch)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(file, destination)
        elif file.suffix == '.log':
            destination = args.output / (str(file.relative_to(args.scratch)) + '.gz')
            destination.parent.mkdir(parents=True, exist_ok=True)
            with file.open('rb') as source, destination.open('wb') as target:
                with gzip.GzipFile(filename='', mode='wb', fileobj=target, mtime=0) as zipped:
                    shutil.copyfileobj(source, zipped)
    hashes = {str(file.relative_to(args.output)): sha(file) for file in sorted(args.output.rglob('*')) if file.is_file() and file.name != 'sha256.json'}
    (args.output / 'sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')
    print('Archived', len(hashes), 'files', flush=True)


if __name__ == '__main__':
    main()
