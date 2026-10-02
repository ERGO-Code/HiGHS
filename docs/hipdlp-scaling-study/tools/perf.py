"""Deterministic per-(instance,config) cost: iterations if solved in any CPU run, else PAR2 penalty."""
import csv, collections, math
S="/tmp/claude-1005/-home-yzhou-Github-hipdlp2/1e818617-21b5-47ff-b840-29674cc27652/scratchpad"
CFGS=["A","L","P","R","C","O","AP","AR"]
def load():
    best=collections.defaultdict(dict); reach=collections.defaultdict(int)
    for f in ("cpu_all.csv","cpu_long.csv"):
        for r in csv.DictReader(open(f"{S}/sweep/{f}")):
            it=int(float(r["iterations"])) if r["iterations"] else 0
            i,c=r["instance"],r["config"]
            reach[i]=max(reach[i],it)
            if r["status"]=="OPTIMAL" and r["kkt_pass"]=="True":
                best[i][c]=min(best[i].get(c,10**12),it)
    # trajectory reruns (A,L,P,AP) add more solved observations
    import glob,re,os
    F=re.compile(r"finished: status ([a-z ]+), (\d+) iterations")
    for p in glob.glob(f"{S}/sweep/traj/*__*.txt"):
        i,c=os.path.basename(p)[:-4].split("__")
        m=F.search(open(p).read())
        if not m: continue
        it=int(m[2]); reach[i]=max(reach[i],it)
        if m[1]=="optimal": best[i][c]=min(best[i].get(c,10**12),it)
    insts=sorted(i for i in best if best[i])
    def cost(i,c): return best[i].get(c, 2*max(reach[i],1))
    def solved(i,c): return c in best[i]
    return insts,cost,solved
def sgm(v,shift=1000.0): return math.exp(sum(math.log(x+shift) for x in v)/len(v))-shift
