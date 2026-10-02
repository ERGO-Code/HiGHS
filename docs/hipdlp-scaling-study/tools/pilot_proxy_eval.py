import csv, math, sys, collections, random, re
sys.path.insert(0,".")
from perf import *
insts,cost,solved=load()
S="/tmp/claude-1005/-home-yzhou-Github-hipdlp2/1e818617-21b5-47ff-b840-29674cc27652/scratchpad"
Kp=collections.defaultdict(dict)   # (inst,pilotcfg,k) -> cfg -> Kw
for r in csv.DictReader(open(f"{S}/diag/proxy_pilot.csv")):
    i,c0,k=r["instance"].split("__")
    Kp[(i,c0,int(k))][r["config"]]={"Kw":float(r["Kw"]),"Kbal":float(r["Kbal"])}
def pilot_solved(i,c0,k):
    t=open(f"{S}/sweep/pilots/{i}__{c0}__{k}.txt.status").read()
    m=re.search(r'"status": "([A-Z_]+)"',t); it=re.search(r'"iterations": (\d+)',t)
    return (m and m[1]=="OPTIMAL" and '"kkt_pass": true' in t), (int(it[1]) if it else k)
def rep(name,pick):
    v=[];s=0
    for i in insts:
        c,cst=pick(i); v.append(cst); s+=solved(i,c)
    print(f"  {name:58s} solved {s:3d}/{len(insts)}  SGM {sgm(v):8.0f}")
print("baselines:")
rep("always A",lambda i:("A",cost(i,"A"))); rep("always L",lambda i:("L",cost(i,"L")))
rep("oracle L,A",lambda i:min((("L",cost(i,"L")),("A",cost(i,"A"))),key=lambda t:t[1]))
for key in ("Kw","Kbal"):
  for c0 in ("L","A"):
    for k in (200,1000,4000):
      for mode in ("continue-if-same","always-restart"):
        def pick(i,c0=c0,k=k,mode=mode,key=key):
            ok,it=pilot_solved(i,c0,k)
            if ok: return c0,it
            kk=Kp.get((i,c0,k))
            if not kk or "L" not in kk or "A" not in kk: return c0,cost(i,c0)
            ch=min(["L","A"],key=lambda c:kk[c][key])
            if ch==c0 and mode=="continue-if-same": return ch,cost(i,ch)
            return ch,k+cost(i,ch)
        rep(f"{key} pilot {c0} k={k:<4d} {mode}",pick)
# accuracy of the pilot estimate on L vs A
print("\naccuracy (faster of L/A, >=1.3x apart):")
for c0 in ("L","A"):
  for k in (200,1000,4000):
    n=ok=0
    for i in insts:
        a,l=cost(i,"A"),cost(i,"L")
        kk=Kp.get((i,c0,k))
        if not kk or max(a,l)<1.3*min(a,l): continue
        n+=1; ok+=(("A" if a<l else "L")==min(["L","A"],key=lambda c:kk[c]["Kw"]))
    print(f"  Kw pilot {c0} k={k}: {ok}/{n}")
