import csv, math, sys, collections, random
sys.path.insert(0,".")
from perf import *
insts,cost,solved=load()
PROXY=sys.argv[1] if len(sys.argv)>1 else "../diag/proxy_true.csv"
K=collections.defaultdict(dict)
for r in csv.DictReader(open(PROXY)):
    K[r["instance"]][r["config"]]={"Kbal":float(r["Kbal"]),"Kw":float(r["Kw"])}
insts=[i for i in insts if i in K and all(c in K[i] for c in CFGS)]
print(len(insts),"instances with proxies for all 8 configs")
def spearman(a,b):
    def rk(v):
        o=sorted(range(len(v)),key=lambda k:v[k]); r=[0]*len(v)
        for pos,k in enumerate(o): r[k]=pos
        return r
    ra,rb=rk(a),rk(b); n=len(a); mu=(n-1)/2
    num=sum((x-mu)*(y-mu) for x,y in zip(ra,rb)); den=sum((x-mu)**2 for x in ra)
    return num/den
for key in ("Kbal","Kw"):
    rho=[spearman([math.log(K[i][c][key]) for c in CFGS],[cost(i,c) for c in CFGS]) for i in insts]
    print(f"{key}: mean within-instance Spearman(proxy, iterations) over 8 configs = {sum(rho)/len(rho):.3f}; median {sorted(rho)[len(rho)//2]:.3f}")
def rep(name,pick):
    v=[cost(i,pick(i)) for i in insts]; s=sum(solved(i,pick(i)) for i in insts)
    print(f"  {name:44s} solved {s:3d}/{len(insts)}  SGM {sgm(v):8.0f}")
print("fixed / oracle:")
for c in ("A","L","P"): rep("always "+c,lambda i,c=c:c)
rep("oracle L,A",lambda i:min(["L","A"],key=lambda c:cost(i,c)))
rep("oracle all 8",lambda i:min(CFGS,key=lambda c:cost(i,c)))
for key in ("Kbal","Kw"):
    print(f"select by argmin {key} (true x*,y*):")
    for cands in (["L","A"],["L","A","P"],["L","A","AP","P"],CFGS):
        rep("among "+",".join(cands),lambda i,cands=cands:min(cands,key=lambda c:K[i][c][key]))
    # accuracy on L vs A when they differ by >=1.3x
    n=ok=0
    for i in insts:
        a,l=cost(i,"A"),cost(i,"L")
        if max(a,l)>=1.3*min(a,l):
            n+=1; ok+= (("A" if a<l else "L")==min(["L","A"],key=lambda c:K[i][c][key]))
    print(f"  picks the faster of L/A (>=1.3x apart): {ok}/{n}")
