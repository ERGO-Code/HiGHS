import os,subprocess,sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,".")
from perf import *
S="/tmp/claude-1005/-home-yzhou-Github-hipdlp2/1e818617-21b5-47ff-b840-29674cc27652/scratchpad"
CONFIGS={"A":["--scaling","auto"],"L":["--scaling","legacy"],"P":["--scaling","legacy","--scaling-mask","4"],"R":["--scaling","legacy","--scaling-mask","1"],
 "C":["--scaling","curtis_reid"],"O":["--scaling","off"],"AP":["--scaling","auto","--scaling-mask","4"],"AR":["--scaling","auto","--scaling-mask","1"]}
insts,cost,solved=load()
def run(i):
    c=min(CONFIGS,key=lambda c:cost(i,c))
    out=f"{S}/sweep/sols/{i}.txt"
    if os.path.exists(out): return
    cmd=[f"{S}/diag/adapter_x_new","--instance",f"/srv/mps_da/{i}.mps.gz","--implementation","rewrite","--backend","cpu",
         "--kkt-tolerance","1e-4","--time-limit","90","--dump-solution",out+".tmp"]+CONFIGS[c]
    p=subprocess.run(cmd,capture_output=True,text=True,timeout=400)
    ok='"kkt_pass": true' in p.stdout
    if ok: os.rename(out+".tmp",out)
    else:
        try: os.remove(out+".tmp")
        except OSError: pass
    print(i,c,"ok" if ok else "FAIL",flush=True)
with ThreadPoolExecutor(14) as ex: list(ex.map(run,insts))
