#!/usr/bin/env python3
"""Serial GPU comparison with matching call-level wall timers and shared KKT."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--external', type=Path, required=True)
    p.add_argument('--names', required=True)
    p.add_argument('--limit', type=float, default=300)
    args = p.parse_args()
    env = dict(os.environ, CUDA_VISIBLE_DEVICES='0', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    records = []
    output = args.root / 'results.json'
    if output.exists():
        records = json.loads(output.read_text())
    for name in args.names.split(','):
        model = args.root / 'models' / (name + '-canonical.mps')
        with model.open('rb') as file:
            model_hash = hashlib.file_digest(file, 'sha256').hexdigest()
        verifier = subprocess.Popen([str(args.root/'verify-server'), str(model)],
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                    text=True, env=env)
        try:
            dimensions = json.loads(verifier.stdout.readline())
            solvers = ['rewrite', 'hpr', 'cupdlpx', 'cupdlpx-no-gm']
            # Both NESO default-geometric failures are already archived in the
            # previous study. Retest accepted configurations with common timers.
            if name.startswith('neso-'):
                solvers.remove('cupdlpx')
            for solver in solvers:
                if any(r['name'] == name and r['solver'] == solver and r.get('selected') for r in records):
                    continue
                tolerances = ['1e-4'] if solver == 'rewrite' else ['1e-4', '1e-6', '1e-8']
                for tolerance in tolerances:
                    out = args.root / 'runs' / name / (solver + '-' + tolerance)
                    out.mkdir(parents=True, exist_ok=True)
                    prefix = out / 'solution'
                    if solver == 'rewrite':
                        cmd = [str(args.root/'rewrite'), str(model), str(prefix), str(args.limit)]
                    elif solver == 'hpr':
                        cmd = [str(args.external/'hpr-build/solve_mps_file'), '-i', str(model),
                               '--presolver', 'none', '--tol', tolerance, '--cr', 'true',
                               '--time-limit', str(args.limit), '--print-debug-info', 'false']
                    else:
                        cmd = [str(args.external/'cupdlpx-current-build/cupdlpx'), '--no_presolve',
                               '--eps_opt', tolerance, '--eps_feas', tolerance, '--time_limit', str(args.limit)]
                        if solver.endswith('no-gm'):
                            cmd += ['--geo_mean_iter', '0']
                        cmd += [str(model), str(out)]
                    record = dict(name=name, solver=solver, native_tolerance=tolerance,
                                  model=str(model), sha256=model_hash, dimensions=dimensions,
                                  command=cmd, time_limit=args.limit, started=time.time())
                    (out/'command.json').write_text(json.dumps(record, indent=2)+'\n')
                    print('START', name, solver, tolerance, flush=True)
                    start = time.monotonic()
                    with (out/'run.log').open('w') as log:
                        try:
                            run = subprocess.run(cmd, env=dict(env, NESO_SOLUTION_PREFIX=str(prefix)),
                                                 stdout=log, stderr=subprocess.STDOUT,
                                                 timeout=args.limit+180)
                            record['exit_code'] = run.returncode
                        except subprocess.TimeoutExpired:
                            record['external_timeout'] = True
                    record['process_wall_sec'] = time.monotonic()-start
                    log = (out/'run.log').read_text()
                    timer = re.findall(r'BENCH_TIMING (\{[^\n]+\})', log)
                    if timer:
                        record.update(json.loads(timer[-1]))
                        record['read_plus_solve_wall_sec'] = record['read_wall_sec'] + record['solver_wall_sec']
                    if (out/'solution.json').exists():
                        native = json.loads((out/'solution.json').read_text())
                        x, y, z = [str(prefix)+suffix for suffix in ('.primal.txt','.dual.txt','.reduced_cost.txt')]
                    elif (out/(model.stem+'_summary.txt')).exists():
                        summary = dict(line.split(': ',1) for line in (out/(model.stem+'_summary.txt')).read_text().splitlines() if ': ' in line)
                        native = dict(status=summary['Termination Reason'], iterations=int(summary['Iterations Count']),
                                      time=float(summary['Runtime (sec)']))
                        x, y, z = [str(out/(model.stem+suffix)) for suffix in ('_primal_solution.txt','_dual_solution.txt','_reduced_cost.txt')]
                    else:
                        native = dict(status='NO_EXPORT')
                        x = y = z = str(out/'missing')
                    record['native'] = native
                    verification = out/'verification.json'
                    verifier.stdin.write('\t'.join((x,y,z,str(verification)))+'\n')
                    verifier.stdin.flush()
                    assert verifier.stdout.readline().strip() == 'DONE'
                    record['verification'] = json.loads(verification.read_text())
                    passed = record['verification']['kkt_pass']
                    terminal = native['status'] in ('TIME_LIMIT','NO_EXPORT')
                    record['selected'] = passed or terminal or tolerance == tolerances[-1]
                    (out/'result.json').write_text(json.dumps(record, indent=2)+'\n')
                    records.append(record)
                    output.write_text(json.dumps(records, indent=2)+'\n')
                    print('DONE', name, solver, tolerance, native['status'],
                          round(record.get('solver_wall_sec', record['process_wall_sec']),3),
                          'KKT', passed, flush=True)
                    if record['selected']:
                        break
        finally:
            verifier.stdin.close()
            verifier.wait()


if __name__ == '__main__':
    main()
