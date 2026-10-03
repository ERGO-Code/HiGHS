#!/usr/bin/env python3
"""Resume sequential, balanced scaling timings on pre-copied scratch LPs.

The driver reads each LP once, warms the backend, then starts a fresh solve for
each configuration. Never run this concurrently with other GPU benchmarks.
"""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import threading
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adapter', type=Path, required=True)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--names', type=Path, required=True)
    parser.add_argument('--backend', choices=['cpu', 'cuda'], default='cuda')
    parser.add_argument('--configs', default='A,L,S,P')
    parser.add_argument('--seconds', type=float, default=180)
    parser.add_argument('--iterations', type=int, default=2147483647)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    names = args.names.read_text().split()
    configs = args.configs.split(',')
    all_results = []
    for index, name in enumerate(names):
        order = configs[index % len(configs):] + configs[:index % len(configs)]
        todo = []
        for config in order:
            result_path = args.output / f'{name}__{config}.json'
            if result_path.exists():
                all_results.append(json.loads(result_path.read_text()))
            else:
                todo.append(config)
        if not todo:
            continue
        prefix = args.output / f'{name}__batch'
        command = [str(args.adapter.resolve()),
                   str(args.scratch / 'lp' / f'{name}.mps'), name,
                   args.backend, ','.join(todo), str(args.seconds),
                   str(args.iterations)]
        Path(str(prefix) + '.command.json').write_text(json.dumps(command) + '\n')
        print(f'BEGIN {index + 1}/{len(names)} {name} {todo}', flush=True)
        started = time.monotonic()
        process = subprocess.Popen(command, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True,
                                   bufsize=1)
        # Covers input loading and verification outside the solver deadline.
        watchdog = threading.Timer(len(todo) * args.seconds + 600, process.kill)
        watchdog.daemon = True
        watchdog.start()
        current = None
        section = []
        model = None
        with Path(str(prefix) + '.log.gz').open('wb') as raw:
            with gzip.GzipFile(fileobj=raw, mode='wb') as log:
                for line in process.stdout:
                    log.write(line.encode())
                    if line.startswith('LARGE_MODEL '):
                        model = json.loads(line[len('LARGE_MODEL '):])
                        (args.output / f'{name}__model.json').write_text(
                            json.dumps(model, indent=2) + '\n')
                    elif line.startswith('LARGE_BEGIN '):
                        current = line.split()[1]
                        section = []
                    elif line.startswith('LARGE_RESULT '):
                        result = json.loads(line[len('LARGE_RESULT '):])
                        result['command'] = command
                        result['guard_rejected'] = any(
                            'Curtis-Reid rejected' in x for x in section)
                        result['selector_log'] = [x.strip() for x in section
                            if 'pilot' in x.lower() or 'selector' in x.lower()]
                        result['model'] = model
                        (args.output / f'{name}__{current}.json').write_text(
                            json.dumps(result, indent=2) + '\n')
                        all_results.append(result)
                        print(f'RESULT {name} {current} {result["status"]} '
                              f'KKT={result["kkt_pass"]} '
                              f'{result["solver_wall_sec"]:.3f}s', flush=True)
                    elif current is not None:
                        section.append(line)
        returncode = process.wait()
        watchdog.cancel()
        (args.output / f'{name}__process.json').write_text(json.dumps({
            'command': command, 'returncode': returncode,
            'process_wall_sec': time.monotonic() - started}, indent=2) + '\n')
        if returncode:
            raise RuntimeError(f'{name}: driver failed ({returncode}); inspect log')
        if len([x for x in all_results if x['name'] == name]) != len(configs):
            raise RuntimeError(f'{name}: missing result')
        (args.output / 'results.json').write_text(
            json.dumps(all_results, indent=2) + '\n')
    print(f'COMPLETE {len(all_results)} solves', flush=True)


if __name__ == '__main__':
    main()
