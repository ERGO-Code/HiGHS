import os,subprocess,sys
from concurrent.futures import ThreadPoolExecutor
S="/tmp/claude-1005/-home-yzhou-Github-hipdlp2/1e818617-21b5-47ff-b840-29674cc27652/scratchpad"
sys.path.insert(0,f"{S}/sweep")
CONFIGS={"A":["--scaling","auto"],"L":["--scaling","legacy"],"P":["--scaling","legacy","--scaling-mask","4"],"AP":["--scaling","auto","--scaling-mask","4"]}
insts=[l.strip() for l in open(f"{S}/sweep/informative.txt") if l.strip()]
os.makedirs(f"{S}/sweep/traj",exist_ok=True)
def run(job):
    i,c=job
    out=f"{S}/sweep/traj/{i}__{c}.txt"
    if os.path.exists(out): return
    cmd=[f"{S}/diag/adapter_x","--instance",f"/srv/mps_da/{i}.mps.gz","--implementation","rewrite","--backend","cpu","--kkt-tolerance","1e-4","--time-limit","30","--log"]+CONFIGS[c]
    p=subprocess.run(cmd,capture_output=True,text=True,timeout=200)
    open(out,"w").write(p.stdout)
jobs=[(i,c) for i in insts for c in CONFIGS]
with ThreadPoolExecutor(12) as ex: list(ex.map(run,jobs))
open(f"{S}/sweep/traj.done","w").write("done")
