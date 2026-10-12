"""Independent check + kill test for the Z16 Construction-A rank-96 candidate.
Lattice: L = {(u, v)/4 : u, v in Z^48, v = u A (mod 16)},  norm(x) = (|u|^2 + |v|^2)/16.
Integer embedding basis rows: [I | A] and [0 | 16 I].
Usage: z16_check_and_kill.py candidate.json runs top_block out.json"""
import sys, json, time, random
import numpy as np
import flint
from fpylll import IntegerMatrix, LLL, BKZ
from fpylll.algorithms.bkz2 import BKZReduction

cand = json.load(open(sys.argv[1]))
runs = int(sys.argv[2]); top = int(sys.argv[3]); out = sys.argv[4]
A = np.array(cand["A"], dtype=np.int64)
n = 48
I = np.eye(n, dtype=np.int64)

# --- structural checks (exact integers) ---
AAt = A @ A.T
ok_sd = ((AAt + I) % 16 == 0).all()
rownorm = 1 + (A * A).sum(1)
ok_t2 = (rownorm % 32 == 0).all()
G = np.block([[(I + AAt) // 16, A], [A.T, 16 * I]])
ok_int = ((I + AAt) % 16 == 0).all()
ok_even = (np.diag(G) % 2 == 0).all()
det = int(flint.fmpz_mat(G.tolist()).det())
ok_pd = np.linalg.eigvalsh(G.astype(float)).min() > 0
print(f"A A^T = -I mod 16: {ok_sd}; rows of [I|A] norm = 0 mod 32: {ok_t2}; Gram integral: {ok_int}; "
      f"even: {ok_even}; det = {det}; positive definite: {ok_pd}", flush=True)
assert ok_sd and ok_t2 and ok_int and ok_even and det == 1 and ok_pd

def check_vec(y):
    u, v = np.array(y[:n]), np.array(y[n:])
    member = ((v - u @ A) % 16 == 0).all()
    return member, (u @ u + v @ v) / 16

wu, wv = cand["short_witness_u"], cand["short_witness_v"]
m, nm = check_vec(list(wu) + list(wv))
print(f"supplied witness: in lattice {m}, norm {nm}", flush=True)

basis = np.vstack([np.hstack([I, A]), np.hstack([np.zeros((n, n), dtype=np.int64), 16 * I])])
results = []
best = None
for r in range(runs):
    rng = random.Random(7000 + r)
    M = IntegerMatrix.from_matrix(basis.tolist())
    for _ in range(4 * 96):
        i, j = rng.sample(range(96), 2)
        M[i].addmul(M[j], rng.choice([-1, 1]))
    LLL.reduction(M)
    t = time.time()
    hist = []
    for bs in [20, 30, 40, 45, 50, 55, 60]:
        if bs > top:
            break
        BKZReduction(M)(BKZ.Param(block_size=bs, max_loops=4, flags=BKZ.AUTO_ABORT | BKZ.MAX_LOOPS))
        rows = [[M[i, j] for j in range(96)] for i in range(96)]
        k = min(range(96), key=lambda i: sum(x * x for x in rows[i]))
        nm = sum(x * x for x in rows[k]) / 16
        hist.append((bs, nm, round(time.time() - t)))
        print(f"run {r} BKZ-{bs}: shortest norm {nm}  ({time.time()-t:.0f}s)", flush=True)
        if best is None or nm < best[0]:
            best = (nm, rows[k])
        if nm < 10:
            break
    results.append(hist)
    if best[0] < 10:
        break
mem, nm = check_vec(best[1])
print(f"BEST: norm {nm}, exact membership {mem}", flush=True)
json.dump(dict(structure_ok=True, runs=results, best_norm=nm, best_vector=[int(x) for x in best[1]],
               best_in_lattice=bool(mem), killed=bool(nm < 10 and mem)), open(out, "w"))
