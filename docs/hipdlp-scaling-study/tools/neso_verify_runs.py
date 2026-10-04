#!/usr/bin/env python3
"""Verify completed exports against the common hipdlp-original-kkt/2 contract."""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--scratch', required=True, type=Path)
    p.add_argument('--verifier', required=True, type=Path)
    p.add_argument('--reported-verifier', type=Path)
    args = p.parse_args()
    records = []
    for process_file in sorted((args.scratch / 'runs').glob('*-canonical/*/process.json')):
        out = process_file.parent
        record = json.loads(process_file.read_text())
        solver = record['solver']
        if solver == 'rewrite-profile':
            continue
        if not any((out / name).exists() for name in
                   ('solution.json', 'hipdlp_summary.json', record['name'] + '_summary.txt')):
            record['verification'] = dict(kkt_pass=False, unverified=True,
                                          reason='No complete native solution export')
            records.append(record)
            print('NO EXPORT', out, flush=True)
            (args.scratch / 'verified-results.json').write_text(json.dumps(records, indent=2) + '\n')
            continue
        if (out / 'solution.json').exists():
            native = json.loads((out / 'solution.json').read_text())
            primal, dual = out / 'solution.primal.txt', out / 'solution.dual.txt'
            costs = out / 'solution.reduced_cost.txt'
        elif solver.startswith('cupdlpx'):
            data = dict(line.split(': ', 1) for line in (out / (record['name'] + '_summary.txt')).read_text().splitlines() if ': ' in line)
            native = dict(status=data['Termination Reason'], iterations=int(data['Iterations Count']), time=float(data['Runtime (sec)']), objective=float(data['Primal Objective Value']))
            primal = out / (record['name'] + '_primal_solution.txt')
            dual = out / (record['name'] + '_dual_solution.txt')
            costs = out / (record['name'] + '_reduced_cost.txt')
        else:
            data = json.loads((out / 'hipdlp_summary.json').read_text())
            native = dict(status=data['status'], iterations=data['iterations'], time=data['solver_wall_sec'], objective=data['objective'])
            lines = (out / 'solution').read_text().splitlines()
            n, m = map(int, lines[0].split())
            assert len(lines) == 1 + 2 * n + m
            primal, dual = out / 'solution.primal.txt', out / 'solution.dual.txt'
            primal.write_text('\n'.join(lines[1:1+n]) + '\n')
            dual.write_text('\n'.join(lines[1+2*n:]) + '\n')
            costs = out / 'solution.reduced_cost.txt'
            costs.write_text('\n'.join(lines[1+n:1+2*n]) + '\n')
        verification_file = out / 'verification.json'
        cmd = [str(args.verifier), '--instance', record['model'], '--primal-values', str(primal), '--row-dual-values', str(dual), '--solver', solver, '--status', native['status'], '--solver-time', str(native['time']), '--iterations', str(native['iterations']), '--tolerance', '1e-4', '--output', str(verification_file)]
        if not verification_file.exists():
            (out / 'verification-command.json').write_text(json.dumps(cmd, indent=2) + '\n')
            with (out / 'verification.log').open('w') as log:
                result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode != 0:
                print('VERIFY ERROR', out, result.returncode, flush=True)
                continue
        verified = json.loads(verification_file.read_text())
        # Metadata does not affect the verifier's numerical calculations.
        verified['solver_time_sec'] = native['time']
        verified['iterations'] = native['iterations']
        verified['solver'] = solver
        verification_file.write_text(json.dumps(verified, indent=2) + '\n')
        record.update(native=native, verification_derived=verified)
        if args.reported_verifier and costs.exists():
            reported_file = out / 'verification-reported.json'
            reported_cmd = [str(args.reported_verifier), record['model'], str(primal), str(dual), str(costs), str(reported_file)]
            if not reported_file.exists():
                (out / 'verification-reported-command.json').write_text(json.dumps(reported_cmd, indent=2) + '\n')
                with (out / 'verification-reported.log').open('w') as log:
                    checked = subprocess.run(reported_cmd, stdout=log, stderr=subprocess.STDOUT)
                if checked.returncode != 0:
                    print('REPORTED VERIFY ERROR', out, checked.returncode, flush=True)
                    continue
            verified = json.loads(reported_file.read_text())
        record['verification'] = verified
        records.append(record)
        print(record['name'], solver, native['status'], native['iterations'], native['time'], 'KKT', verified['kkt_pass'], flush=True)
        (args.scratch / 'verified-results.json').write_text(json.dumps(records, indent=2) + '\n')


if __name__ == '__main__':
    main()
