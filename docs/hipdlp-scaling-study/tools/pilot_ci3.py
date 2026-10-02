import sys,math,random
sys.argv=[sys.argv[0]]
exec(open("pilot_proxy_eval.py").read().split('print("baselines:")')[0])
random.seed(7)
def gm(x): return math.exp(sum(math.log(v) for v in x)/len(x))
def two_pilot(k,how,first="A"):
    other="L" if first=="A" else "A"
    def pick(i):
        # run `first` for k; if it already solved, done. Otherwise run `other` for k too.
        ok,it=pilot_solved(i,first,k)
        if ok: return first,it
        ok2,it2=pilot_solved(i,other,k)
        if ok2: return other,k+it2
        a=Kp[(i,"A",k)]; l=Kp[(i,"L",k)]
        if how=="cross": score={"A":l["A"]["Kw"],"L":a["L"]["Kw"]}
        elif how=="geo": score={c:math.sqrt(a[c]["Kw"]*l[c]["Kw"]) for c in ("A","L")}
        elif how=="self": score={"A":a["A"]["Kw"],"L":l["L"]["Kw"]}
        ch=min(("A","L"),key=lambda c:score[c])
        # both pilots ran k iterations; the chosen one continues, the other is wasted
        return ch,cost(i,ch)+k
    return pick
rows=[]
for k in (200,1000,4000):
  for how in ("cross","geo","self"):
    pick=two_pilot(k,how)
    rA=[];rL=[];worst=[];s=0
    for i in insts:
        c,cs=pick(i); s+=solved(i,c)
        rA.append((cs+1000)/(cost(i,"A")+1000)); rL.append((cs+1000)/(cost(i,"L")+1000)); worst.append((cs/cost(i,"A"),i,c))
    b=sorted(gm([rA[random.randrange(len(rA))] for _ in rA]) for _ in range(3000))
    worst.sort(reverse=True)
    print(f"two pilots k={k:<4d} {how:5s}: solved {s}/{len(insts)}  vs A GM {gm(rA):.3f} CI[{b[75]:.3f},{b[2925]:.3f}]  vs L GM {gm(rL):.3f}  >=3x worse than A: {sum(w[0]>=3 for w in worst)}  worst: "+", ".join(f"{w[1]}({w[2]}) {w[0]:.1f}x" for w in worst[:2]))
