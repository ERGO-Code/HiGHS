#!/usr/bin/env python3
"""Run NESO comparisons serially on one GPU, retaining commands and full logs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--models', type=Path, required=True)
    p.add_argument('--solvers', default='cupdlpx,hpr-cr,hpr-no-cr,rewrite-auto,rewrite-legacy')
    p.add_argument('--names', default='neso-2245,neso-2005')
    p.add_argument('--time-limit', type=float, default=600)
    p.add_argument('--hpr-tolerance', default='1e-4')
    p.add_argument('--cup-tolerance', default='1e-4')
    p.add_argument('--wait-for-pid', type=int)
    args = p.parse_args()
    if args.wait_for_pid:
        print('WAIT FOR GPU RUNNER', args.wait_for_pid, flush=True)
        while True:
            try:
                os.kill(args.wait_for_pid, 0)
            except ProcessLookupError:
                break
            time.sleep(5)
    env = dict(os.environ, CUDA_VISIBLE_DEVICES='0', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for name in args.names.split(','):
        model = args.models / (name + '.mps')
        sha = hashlib.sha256(model.read_bytes()).hexdigest()
        for solver in args.solvers.split(','):
            out = args.scratch / 'runs' / name / solver
            out.mkdir(parents=True, exist_ok=True)
            if solver == 'rewrite-profile':
                cmd = [str(args.scratch / 'rewrite-profile'), str(model), 'legacy']
            elif solver.startswith('cupdlpx'):
                extra = ['--geo_mean_iter', '0'] if 'no-gm' in solver else []
                cmd = [str(args.scratch / 'cupdlpx-current-build/cupdlpx'), '--no_presolve', '--eps_opt', args.cup_tolerance, '--eps_feas', args.cup_tolerance, '--time_limit', str(args.time_limit), '--verbose'] + extra + [str(model), str(out)]
            elif solver.startswith('julia-'):
                cmd = ['julia', '--project=' + str(args.scratch / 'HPR-LP-current'), str(Path(__file__).with_name('neso_hpr_julia.jl')), str(model), str(out / 'solution'), str(solver.startswith('julia-cr')).lower(), args.hpr_tolerance, str(args.time_limit)]
            elif solver.startswith('hpr-'):
                cmd = [str(args.scratch / 'hpr-build/solve_mps_file'), '-i', str(model), '--presolver', 'none', '--tol', args.hpr_tolerance, '--time-limit', str(args.time_limit), '--cr', str(solver.startswith('hpr-cr')).lower(), '--print-debug-info', 'true']
            else:
                cmd = [str(args.scratch / 'rewrite-adapter'), '--instance', str(model), '--implementation', 'rewrite', '--backend', 'cuda', '--step-size', 'pid', '--presolve', 'off', '--scaling', solver.removeprefix('rewrite-'), '--kkt-tolerance', '1e-4', '--time-limit', str(args.time_limit), '--iteration-limit', '2147483647', '--output-dir', str(out), '--dump-solution', str(out / 'solution'), '--log']
            record = dict(name=name, solver=solver, model=str(model), sha256=sha, command=cmd, time_limit=args.time_limit, started=time.time(), hpr_tolerance=args.hpr_tolerance, cup_tolerance=args.cup_tolerance)
            (out / 'command.json').write_text(json.dumps(record, indent=2) + '\n')
            print('START', name, solver, flush=True)
            start = time.monotonic()
            with (out / 'run.log').open('w') as log:
                try:
                    run_env = dict(env, NESO_SOLUTION_PREFIX=str(out / 'solution'))
                    result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=run_env, timeout=args.time_limit + 300)
                    record['exit_code'] = result.returncode
                except subprocess.TimeoutExpired:
                    record['external_timeout'] = True
                except OSError as error:
                    record['launch_error'] = str(error)
            record['process_wall_sec'] = time.monotonic() - start
            (out / 'process.json').write_text(json.dumps(record, indent=2) + '\n')
            print('DONE', name, solver, record.get('exit_code'), round(record['process_wall_sec'], 2), flush=True)


if __name__ == '__main__':
    main()
