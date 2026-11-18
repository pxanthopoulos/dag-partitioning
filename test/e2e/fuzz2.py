#!/usr/bin/env python3
# Like fuzz.py, on DAGs with parallel chains that join late, with big and small
# tensors, so that the order of the partitions has freedom and the limit binds.
# Usage: fuzz2.py <dir of dag-test> <work dir> <graphs> [seed]
import os, re, subprocess, sys, random
os.environ["DAG_PARTITIONING_RANDOM_SEED"] = "42"
B = os.path.abspath(sys.argv[1]) + "/"
out = sys.argv[2]; os.makedirs(out, exist_ok=True); os.chdir(out)
N = int(sys.argv[3]); random.seed(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
def gen(path):
    chains = random.randint(2, 4); nodes = []; edges = []; w = {}; idx = 0; tails = []
    for c in range(chains):
        prev = None
        for k in range(random.randint(2, 4)):
            w[idx] = random.choice([16, 32, 64, 256, 1024, 1024, 4096]); 
            if prev is not None: edges.append((prev, idx))
            prev = idx; idx += 1
        tails.append(prev)
    prevj = None
    random.shuffle(tails)
    for t in tails:
        j = idx; idx += 1; w[j] = random.choice([16, 64, 256, 1024])
        edges.append((t, j))
        if prevj is not None: edges.append((prevj, j))
        if random.random() < 0.3 and len(tails) > 1:
            o = random.choice(tails)
            if o != t: edges.append((o, j))
        prevj = j
    edges = sorted(set(edges))
    with open(path, "w") as f:
        f.write("strict digraph  {\n")
        for n in range(idx): f.write(f'{n} [weight={w[n]}];\n')
        for a, b in edges: f.write(f'{a} -> {b} [weight={w[a]}];\n')
        f.write("}\n")
    return idx
def run(P, dot, *lim):
    r = subprocess.run([B+"dag-test", str(P), dot, "FORB", "GGG", "BOUNDARYFM", "0", "100", "10", "1", "transfer-memory", *map(str, lim)], capture_output=True, text=True, timeout=900)
    m = re.search(r"cost (\d+) peak (\d+)", r.stderr)
    return r.returncode, (int(m.group(1)), int(m.group(2))) if m else None, r.stderr[-500:]
tot = fail = bind = cost = 0; ex = []
for i in range(N):
    dot = f"h{i}.dot"; n = gen(dot); P = random.randint(3, min(8, n - 1))
    rc, free, e = run(P, dot)
    if rc or not free: print("FAIL none", i, n, P, e); fail += 1; continue
    rc, low, e = run(P, dot, 1)
    if rc or not low: print("FAIL limit1", i, n, P, e); fail += 1; continue
    tot += 2
    if low[1] > free[1]: print("LOWEST ABOVE FREE", i, free, low); fail += 1
    for L in sorted({low[1], (low[1]+free[1])//2, free[1]-1, free[1]}):
        if L < low[1]: continue
        rc, res, e = run(P, dot, L); tot += 1
        if rc or not res: print("FAIL", i, n, P, L, e); fail += 1; continue
        if res[1] > L: print("OVER LIMIT", i, P, L, res); fail += 1
        if res[0] < free[0]: print("BEATS UNCONSTRAINED", i, P, L, res, free); fail += 1
    if low[1] < free[1]:
        bind += 1; cost += low[0] > free[0]
        if low[0] > free[0] and len(ex) < 3: ex.append((i, P, "free (moved,peak)", free, "lowest", low))
print(f"{tot} solves on {N} graphs, {fail} failures; limit binds on {bind} graphs; meeting it costs extra transfers on {cost}")
for x in ex: print(" ", x)
sys.exit(1 if fail else 0)
