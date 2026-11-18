#!/usr/bin/env python3
# CP-SAT against brute force (dag-test verifies internally and fails on any
# mismatch) on random DAGs, with no limit, a 1-byte limit and limits in between.
# Usage: fuzz.py <dir of dag-test and rand-dag> <work dir> <graphs> [seed]
import os, re, subprocess, sys, random
os.environ["DAG_PARTITIONING_RANDOM_SEED"] = "42"
B = os.path.abspath(sys.argv[1]) + "/"
out = sys.argv[2]; os.makedirs(out, exist_ok=True); os.chdir(out)
N = int(sys.argv[3]); random.seed(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
def run(P, dot, *lim):
    r = subprocess.run([B+"dag-test", str(P), dot, "FORB", "GGG", "BOUNDARYFM", "0", "100", "10", "1", "transfer-memory", *map(str, lim)], capture_output=True, text=True, timeout=900)
    m = re.search(r"cost (\d+) peak (\d+)", r.stderr)
    return r.returncode, (int(m.group(1)), int(m.group(2))) if m else None, r.stderr[-600:]
tot = fail = bind = lowered = 0
for i in range(N):
    size, ratio, P = random.choice([20, 30, 50, 80]), random.choice([10, 20, 40]), random.choice([3, 4, 5, 6, 7, 8])
    dot = f"g{i}.dot"
    subprocess.run([B+"rand-dag", str(size), str(ratio), "0", dot], capture_output=True, timeout=120)
    rc, free, e = run(P, dot)
    if rc or not free: print("FAIL none", i, size, ratio, P, e); fail += 1; continue
    rc, low, e = run(P, dot, 1)
    if rc or not low: print("FAIL limit1", i, size, ratio, P, e); fail += 1; continue
    tot += 2
    assert low[1] <= free[1], (free, low)
    lims = sorted({low[1], (low[1]+free[1])//2, free[1]-1, free[1]} - {0})
    for L in lims:
        if L < low[1]: continue
        rc, res, e = run(P, dot, L); tot += 1
        if rc or not res: print("FAIL", i, size, ratio, P, L, e); fail += 1; continue
        if res[1] > L: print("OVER", i, P, L, res); fail += 1
        if res[0] < free[0]: print("CHEAPER THAN UNCONSTRAINED", i, P, L, res, free); fail += 1
    if low[1] < free[1]:
        bind += 1; lowered += low[0] > free[0]
print(f"{tot} solves on {N} graphs, {fail} failures; on {bind} graphs the limit binds, and in {lowered} of them meeting it costs extra transfers")
sys.exit(1 if fail else 0)
