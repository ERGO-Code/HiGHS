#!/usr/bin/env python3
"""One fixed-1e-4 run per model and solver; presolve ON, Ruiz10+PC only."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time


def sha256(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--external', type=Path, required=True)
    p.add_argument('--names', default='')
    p.add_argument('--limit', type=float, default=300)
    args = p.parse_args()
    env = dict(os.environ, CUDA_VISIBLE_DEVICES='0', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    output = args.root/'results.json'
    records = json.loads(output.read_text()) if output.exists() else []
    names = args.names.split(',') if args.names else json.loads((args.root/'names.json').read_text())
    for name in names:
        if all(any(r['name']==name and r['solver']==s for r in records) for s in ['legacy','rewrite','hpr','cupdlpx']):
            continue
        model = args.root/'models'/(name+'-canonical.mps')
        model_hash = sha256(model)
        verifier = subprocess.Popen([str(args.root/'verify-server'), str(model)],
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                    text=True, env=env)
        try:
            dimensions = json.loads(verifier.stdout.readline())
            for solver in ['legacy','rewrite','hpr','cupdlpx']:
                if any(r['name']==name and r['solver']==solver for r in records):
                    continue
                out = args.root/'runs'/name/solver
                out.mkdir(parents=True, exist_ok=True)
                prefix = out/'solution'
                if solver in ['legacy','rewrite']:
                    cmd = [str(args.root/'hipdlp-presolve'),str(model),solver,str(prefix),str(args.limit)]
                elif solver == 'hpr':
                    cmd = [str(args.external/'hpr-build/solve_mps_file'),'-i',str(model),
                           '--presolver','pslp','--tol','1e-4','--cr','false',
                           '--ruiz','true','--pock','true','--bc','false',
                           '--time-limit',str(args.limit),'--print-debug-info','false']
                else:
                    cmd = [str(args.external/'cupdlpx-current-build/cupdlpx'),
                           '--eps_opt','1e-4','--eps_feas','1e-4','--time_limit',str(args.limit),
                           '--geo_mean_iter','0','--l_inf_ruiz_iter','10',
                           '--pock_chambolle_alpha','1','--no_bound_obj_rescaling',str(model),str(out)]
                record = dict(name=name,solver=solver,native_tolerance=1e-4,scaling='Ruiz10+PC(alpha=1)',
                              presolve=True,model=str(model),sha256=model_hash,dimensions=dimensions,
                              executable_sha256=sha256(cmd[0]),command=cmd,time_limit=args.limit,
                              external_wall_limit=2*args.limit+180,started=time.time())
                (out/'command.json').write_text(json.dumps(record,indent=2)+'\n')
                print('START', name, solver, flush=True)
                start = time.monotonic()
                with (out/'run.log').open('w') as log:
                    try:
                        run = subprocess.run(cmd,env=dict(env,NESO_SOLUTION_PREFIX=str(prefix)),
                                             stdout=log,stderr=subprocess.STDOUT,
                                             timeout=record['external_wall_limit'])
                        record['exit_code'] = run.returncode
                    except subprocess.TimeoutExpired:
                        record['external_timeout'] = True
                record['process_wall_sec'] = time.monotonic()-start
                log = (out/'run.log').read_text()
                timer = re.findall(r'BENCH_TIMING (\{[^\n]+\})',log)
                if timer:
                    record.update(json.loads(timer[-1]))
                    record['read_plus_solve_wall_sec'] = record['read_wall_sec']+record['solver_wall_sec']
                if (out/'solution.json').exists():
                    native = json.loads((out/'solution.json').read_text())
                    x,y,z = [str(prefix)+s for s in ['.primal.txt','.dual.txt','.reduced_cost.txt']]
                elif (out/(model.stem+'_summary.txt')).exists():
                    summary = dict(line.split(': ',1) for line in (out/(model.stem+'_summary.txt')).read_text().splitlines() if ': ' in line)
                    native = dict(status=summary['Termination Reason'],iterations=int(summary['Iterations Count']),time=float(summary['Runtime (sec)']))
                    x,y,z = [str(out/(model.stem+s)) for s in ['_primal_solution.txt','_dual_solution.txt','_reduced_cost.txt']]
                else:
                    native = dict(status='EXTERNAL_TIMEOUT' if record.get('external_timeout') else 'NO_EXPORT')
                    x=y=z=str(out/'missing')
                record['native'] = native
                verification = out/'verification.json'
                verifier.stdin.write('\t'.join([x,y,z,str(verification)])+'\n')
                verifier.stdin.flush()
                assert verifier.stdout.readline().strip() == 'DONE'
                record['verification'] = json.loads(verification.read_text())
                record['exports'] = {k:dict(path=v,sha256=sha256(v),bytes=Path(v).stat().st_size) for k,v in zip(['primal','dual','reduced_cost'],[x,y,z]) if Path(v).exists()}
                (out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
                records.append(record)
                temporary = output.with_suffix('.tmp')
                temporary.write_text(json.dumps(records,indent=2)+'\n')
                temporary.replace(output)
                print('DONE',name,solver,native['status'],round(record.get('solver_wall_sec',record['process_wall_sec']),3),
                      'KKT',record['verification']['kkt_pass'],flush=True)
        finally:
            verifier.stdin.close()
            verifier.wait()


if __name__ == '__main__':
    main()
