#!/usr/bin/env python3
"""Copy read-only benchmarks and decode complete gzip members into scratch.

Some benchmark .gz files contain a valid MPS gzip member followed by plain-text
SCIP logs. Python gzip's concatenated-member reader raises on that suffix and
can leave a partially written MPS. Decode and validate each actual gzip member,
ignore only non-gzip trailing text, and record its size and the payload hash.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zlib


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def decode(source, destination):
    temporary = destination.with_suffix('.mps.complete')
    with source.open('rb') as stream, temporary.open('wb') as output:
        pending = stream.read(65536)
        if not pending.startswith(b'\x1f\x8b'):
            output.write(pending)
            shutil.copyfileobj(stream, output)
            members = trailing = 0
        else:
            members = trailing = 0
            while pending:
                if not pending.startswith(b'\x1f\x8b'):
                    trailing = len(pending) + sum(len(b) for b in iter(
                        lambda: stream.read(65536), b''))
                    break
                decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
                members += 1
                while True:
                    output.write(decoder.decompress(pending))
                    if decoder.eof:
                        pending = decoder.unused_data
                        if not pending:
                            pending = stream.read(65536)
                        elif len(pending) == 1:
                            pending += stream.read(65536)
                        break
                    pending = stream.read(65536)
                    if not pending:
                        raise ValueError(f'{source}: incomplete gzip member')
    temporary.replace(destination)
    return members, trailing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('/srv/mps_da'))
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--names', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    args = parser.parse_args()
    lp = args.scratch / 'lp'
    lp.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name in args.names.read_text().split():
        source = args.source / f'{name}.mps.gz'
        copy = lp / source.name
        destination = lp / f'{name}.mps'
        previous = digest(destination) if destination.exists() else None
        shutil.copyfile(source, copy)
        members, trailing = decode(copy, destination)
        row = {'name': name, 'source': str(source), 'sha256': digest(copy),
               'mps_sha256': digest(destination), 'gzip_members': members,
               'non_gzip_trailing_bytes': trailing, 'previous_mps_sha256': previous}
        row['previous_copy_identical'] = previous == row['mps_sha256']
        manifest.append(row)
        print(name, 'identical' if row['previous_copy_identical'] else 'CORRECTED',
              'trailing bytes', trailing, flush=True)
        args.manifest.write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
