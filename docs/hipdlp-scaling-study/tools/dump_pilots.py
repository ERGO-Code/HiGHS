import os,subprocess,sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,".")
from perf import *
S="/tmp/claude-1005/-home-yzhou-Github-hipdlp2/1e818617-21b5-47ff-b840-29674cc27652/scratchpad"
CONFIGS={"A":["--scaling","auto"],"L":["--scaling","legacy"]}
insts,cost,solved=load()
jobs=[(i,c,k) for i in insts for c in ("L","A") for k in (200,1000,4000)]
def run(job):
    i,c,k=job
    out=f"{S}/sweep/pilots/{i}__{c}__{k}.txt"
    if os.path.exists(out): return
    cmd=[f"{S}/diag/adapter_x_new","--instance",f"/srv/mps_da/{i}.mps.gz","--implementation","rewrite","--backend","cpu",
         "--kkt-tolerance","1e-4","--time-limit","120","--iteration-limit",str(k),"--dump-solution",out+".tmp"]+CONFIGS[c]
    p=subprocess.run(cmd,capture_output=True,text=True,timeout=600)
    if os.path.exists(out+".tmp"): os.rename(out+".tmp",out)
    with open(out+".status","w") as f: f.write(p.stdout)
with ThreadPoolExecutor(14) as ex: list(ex.map(run,jobs))
print("done")
