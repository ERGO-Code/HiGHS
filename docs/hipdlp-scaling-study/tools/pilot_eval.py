#!/usr/bin/env python3
"""Offline evaluation of 'pilot then pick' scaling selectors from the logged trajectories."""
import glob, math, os, re, sys, collections

S = "/tmp/claude-1005/-home-yzhou-Github-hipdlp2/1e818617-21b5-47ff-b840-29674cc27652/scratchpad"
B = re.compile(r"boundary: iteration (\d+), relative primal ([\d.eE+-]+), dual ([\d.eE+-]+), gap ([\d.eE+-]+)")
F = re.compile(r"finished: status ([a-z ]+), (\d+) iterations, ([\d.]+)s")
EPS = 1e-4
CANDS = sys.argv[1].split(",") if len(sys.argv) > 1 else ["A", "L"]

traj = {}
for p in glob.glob(f"{S}/sweep/traj/*__*.txt"):
    inst, cfg = os.path.basename(p)[:-4].split("__")
    txt = open(p).read()
    b = [(int(m[1]), float(m[2]), float(m[3]), float(m[4])) for m in B.finditer(txt)]
    f = F.search(txt)
    if not f: continue
    traj[(inst, cfg)] = dict(bounds=b, status=f[1], iters=int(f[2]), time=float(f[3]))

insts = sorted({i for i, _ in traj if all((i, c) in traj for c in ["A", "L", "P", "AP"])})
solved = lambda i, c: traj[(i, c)]["status"] == "optimal"
informative = [i for i in insts if any(solved(i, c) for c in ["A", "L", "P", "AP"])]
print(f"{len(insts)} instances with trajectories, {len(informative)} informative")

# per-instance budget (iterations reached by the longest run) for the PAR2 penalty
def budget(i): return max(traj[(i, c)]["iters"] for c in ["A", "L", "P", "AP"])

def cost_of(i, c):
    return traj[(i, c)]["iters"] if solved(i, c) else 2 * budget(i)

def sgm(v, shift=1000.0):
    return math.exp(sum(math.log(x + shift) for x in v) / len(v)) - shift

def E_at(i, c, k):
    """merit at the first boundary >= k (None if run ended earlier)"""
    for it, rp, rd, g in traj[(i, c)]["bounds"]:
        if it >= k: return max(rp, rd, g), it
    return None, None

def merit_log(i, c, k):
    for it, rp, rd, g in traj[(i, c)]["bounds"]:
        if it >= k: return sum(math.log10(max(x, 1e-16)) for x in (rp, rd, g)) / 3.0
    return None

def report(name, pick, overhead_k=0):
    costs, nsolved = [], 0
    for i in informative:
        c, ov = pick(i)
        cst = cost_of(i, c) + ov
        # if the chosen run was unsolved the cost stays the penalty
        costs.append(cst); nsolved += solved(i, c)
    print(f"{name:46s} solved {nsolved:3d}/{len(informative)}  PAR2-SGM(iters) {sgm(costs):10.0f}")

print("\n--- fixed")
for c in ["A", "L", "P", "AP"]:
    report(f"always {c}", lambda i, c=c: (c, 0))
report("ORACLE over A,L", lambda i: (min(["A", "L"], key=lambda c: cost_of(i, c)), 0))
report("ORACLE over A,L,P,AP", lambda i: (min(["A", "L", "P", "AP"], key=lambda c: cost_of(i, c)), 0))

print(f"\n--- pilot selectors over {CANDS}: choose by merit at pilot length k; cost = chosen run + (n-1)*k pilot iterations")
for kind in ("max", "logmean"):
    for k in (1600, 3200, 6400, 12800):
        def pick(i, k=k, kind=kind):
            vals = {}
            for c in CANDS:
                # a candidate that already solved inside the pilot window ends the pilot
                if solved(i, c) and traj[(i, c)]["iters"] <= k:
                    vals[c] = -1e9
                    continue
                v = E_at(i, c, k)[0] if kind == "max" else merit_log(i, c, k)
                if v is None: v = 1e9
                vals[c] = v
            best = min(CANDS, key=lambda c: vals[c])
            ov = (len(CANDS) - 1) * min(k, traj[(i, best)]["iters"])
            if vals[best] == -1e9: ov = 0
            return best, ov
        report(f"pilot {kind:7s} k={k}", pick)
