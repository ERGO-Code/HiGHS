import hashlib,json,subprocess
from pathlib import Path
root=Path('/tmp/hipdlp-four20'); repo=Path('/home/yzhou/Documents/Github/hipdlp2')
def cmd(*c,cwd=repo): return subprocess.check_output(c,cwd=cwd,text=True).strip()
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
extern=Path('/tmp/neso-gpu-comparison')
meta=dict(date=cmd('date','--iso-8601=seconds'),hipdlp_commit=cmd('git','rev-parse','HEAD'),
          hpr_commit=cmd('git','rev-parse','HEAD',cwd=extern/'HPR-LP-C'),
          cupdlpx_commit=cmd('git','rev-parse','HEAD',cwd=extern/'cuPDLPx-current'),
          gpu=cmd('nvidia-smi','--query-gpu=name,driver_version,memory.total','--format=csv,noheader'),
          cpu=cmd('lscpu'),cuda=cmd('/usr/local/cuda-12.8/bin/nvcc','--version'),
          environment=dict(CUDA_VISIBLE_DEVICES='0',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),
          solvers_order=['legacy','rewrite','hpr','cupdlpx'],native_tolerance=1e-4,
          shared_kkt_tolerance=1e-4,native_limit_seconds=300,external_limit_seconds=780,
          scaling=dict(ruiz_iterations=10,pc_alpha=1,cr=False,geometric=False,bound_cost=False,l2=False),
          library_sha256=sha(root/'build/lib/libhighs.so'))
meta['artifacts']={str(p):sha(p) for p in [root/'hipdlp-presolve',root/'verify-server',extern/'hpr-build/solve_mps_file',extern/'cupdlpx-current-build/cupdlpx',repo/'docs/hipdlp-scaling-study/tools/four_solver_presolve.cpp',repo/'docs/hipdlp-scaling-study/tools/four_solver_run.py',repo/'docs/hipdlp-scaling-study/tools/gpu_comparison_verify.cpp']}
(root/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
for source,file in [(extern/'HPR-LP-C','hpr-passive-cli.patch'),(extern/'cuPDLPx-current','cupdlpx-passive-cli.patch')]:
    (root/file).write_text(cmd('git','diff',cwd=source)+'\n')
