#!/usr/bin/env python3
"""Offline evaluation of scaling-selection policies from the sweep CSV + per-config features."""
import csv, math, sys, collections

S = "/tmp/claude-1005/-home-yzhou-Github-hipdlp2/1e818617-21b5-47ff-b840-29674cc27652/scratchpad"
TL = float(sys.argv[2]) if len(sys.argv) > 2 else 30.0
CSV = sys.argv[1] if len(sys.argv) > 1 else f"{S}/sweep/cpu_all.csv"
SHIFT = 1.0
CFGS = ["A", "L", "P", "R", "C", "O", "AP", "AR"]

feat = {}
for r in csv.DictReader(open(f"{S}/diag/feat_all.csv")):
    feat[(r["instance"], r["config"])] = {k: (float(v) if k not in ("instance", "config") else v) for k, v in r.items()}

res = collections.defaultdict(dict)   # inst -> cfg -> (solved, time, iters)
for r in csv.DictReader(open(CSV)):
    ok = r["status"] == "OPTIMAL" and r["kkt_pass"] == "True"
    t = float(r["solver_wall_sec"]) if r["solver_wall_sec"] else TL
    res[r["instance"]][r["config"]] = (ok, t if ok else 2 * TL, int(float(r["iterations"])) if r["iterations"] else 0)

inst = [i for i in sorted(res) if all(c in res[i] for c in CFGS)]
# informative: at least one of the 8 solves it
info = [i for i in inst if any(res[i][c][0] for c in CFGS)]
print(f"{len(inst)} complete instances, {len(info)} solved by at least one config")

def lg(x): return math.log10(x) if x > 0 else 0.0

def metrics(i):
    a, l = feat[(i, "A")], feat[(i, "L")]
    p = feat.get((i, "P"))
    m = {}
    m["amp"] = lg(max(a["ampR_max"], a["ampC_max"]))
    m["ampC"] = lg(a["ampC_max"]); m["ampR"] = lg(a["ampR_max"])
    m["amp_ratio"] = m["amp"] - lg(max(l["ampR_max"], l["ampC_max"]))
    m["span99"] = (a["Rp99_l"] - a["Rp1_l"]) + (a["Cp99_l"] - a["Cp1_l"])
    m["span_full"] = (a["Rmax_l"] - a["Rmin_l"]) + (a["Cmax_l"] - a["Cmin_l"])
    ratios = []
    for k in ("cn", "bn", "bndn"):
        if l[k] > 0 and a[k] > 0: ratios.append(lg(a[k] / l[k]))
    m["norm_ratio"] = max(ratios) if ratios else 0.0
    return m

M = {i: metrics(i) for i in inst}

def sgm(vals):
    return math.exp(sum(math.log(v + SHIFT) for v in vals) / len(vals)) - SHIFT

def evaluate(policy, name):
    ts, solved, worst = [], 0, 0.0
    for i in info:
        c = policy(i)
        ok, t, it = res[i][c]
        ts.append(t); solved += ok
        best = min(res[i][x][1] for x in CFGS)
        worst = max(worst, t / max(best, 0.05))
    print(f"{name:34s} solved {solved:3d}/{len(info)}  PAR2-SGM {sgm(ts):7.3f}  worst-vs-best {worst:8.1f}x")
    return sgm(ts), solved

print("\n--- fixed policies")
for c in CFGS:
    evaluate(lambda i, c=c: c, f"always {c}")
evaluate(lambda i: min(CFGS, key=lambda c: res[i][c][1]), "ORACLE (best of 8)")
evaluate(lambda i: min(["A", "L"], key=lambda c: res[i][c][1]), "ORACLE (best of A,L)")
evaluate(lambda i: min(["A", "L", "P"], key=lambda c: res[i][c][1]), "ORACLE (best of A,L,P)")

print("\n--- guard policies: use A unless metric > T, then fall back")
for fb in ("L", "P"):
    for key in ("amp", "ampC", "ampR", "amp_ratio", "span99", "span_full", "norm_ratio"):
        vals = sorted({round(M[i][key], 1) for i in info})
        best = None
        for T in vals:
            def pol(i, T=T, key=key, fb=fb): return fb if M[i][key] > T else "A"
            ts = [res[i][pol(i)][1] for i in info]
            sol = sum(res[i][pol(i)][0] for i in info)
            v = sgm(ts)
            if best is None or (sol, -v) > (best[1], -best[0]): best = (v, sol, T)
        T = best[2]
        # also show the neso-like threshold stability: value range with solved >= best-1
        print(f"fallback {fb}  metric {key:11s} best T={T:6.1f}  solved {best[1]}  PAR2-SGM {best[0]:.3f}")
