"""Heavier kill test for survivors: several randomized bases, LLL -> BKZ 20/30/40 (and optionally 45).
Reports shortest norm found per run."""
import sys, json, time, random
import numpy as np
sys.path.insert(0, sys.argv[1])
from fpylll import IntegerMatrix, LLL, BKZ
from fpylll.algorithms.bkz2 import BKZReduction
import importlib.util
spec = importlib.util.spec_from_file_location("t96", sys.argv[1] + "/tensor96_lib.py")
t96 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t96)

D = sys.argv[1]
cands = json.load(open(sys.argv[2]))
idxs = [int(x) for x in sys.argv[3].split(",")]
runs = int(sys.argv[4]); top_block = int(sys.argv[5]); out = sys.argv[6]
leech = json.load(open(f"{D}/leech23.json"))
res = []
for pi in idxs:
    P = cands[pi]
    st = [s for s in leech if s["sign"] == P["sign"]][0]
    Y, Gm, det, even, Yamb, p = t96.build(P, st)
    assert det == 1 and even
    scale = 8 * p
    best_overall = None
    for r in range(runs):
        rng = random.Random(1000 * pi + r)
        A = IntegerMatrix.from_matrix([[int(v) for v in row] for row in Yamb])
        n = A.nrows
        # randomize: a few hundred random row operations
        for _ in range(4 * n):
            i, j = rng.sample(range(n), 2)
            c = rng.choice([-1, 1])
            A[i].addmul(A[j], c)
        LLL.reduction(A)
        t = time.time()
        for bs in [20, 30, 40, 45, 50]:
            if bs > top_block:
                break
            BKZReduction(A)(BKZ.Param(block_size=bs, max_loops=4, flags=BKZ.AUTO_ABORT | BKZ.MAX_LOOPS))
            m = min(sum(A[i, j] ** 2 for j in range(A.ncols)) for i in range(n)) / scale
            if m < 10:
                break
        vec = min((sum(A[i, j] ** 2 for j in range(A.ncols)), i) for i in range(n))
        m = vec[0] / scale
        w = [A[vec[1], j] for j in range(A.ncols)]
        print(f"P#{pi} run {r}: shortest norm {m} ({time.time()-t:.0f}s)", flush=True)
        if best_overall is None or m < best_overall[0]:
            best_overall = (m, w)
        if m < 10:
            break
    res.append(dict(P_index=pi, p=P["p"], sign=P["sign"], best=best_overall[0], witness=best_overall[1]))
    print(f"==> P#{pi}: best {best_overall[0]} {'KILLED' if best_overall[0] < 10 else 'survives'}", flush=True)
    json.dump(res, open(out, "w"))
