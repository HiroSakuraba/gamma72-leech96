"""Rank-4 Hermitian unimodular lattices over O = Z[beta], beta = (1+sqrt(-11))/2 (beta^2 - beta + 3 = 0),
from codes over split primes p:  P_C = { x in O^4 : x mod P1 in C1, x mod P2 in C2 = C1^perp },
form h = (1/p) sum x_i conj(y_i).  Goal: Hermitian min >= 3  (trace min 6, extremal 11-modular in dim 8).
Exhaustive over all C1 for small p, random for larger p."""
import sys, json, itertools, random, math
import numpy as np
import flint
sys.path.insert(0, sys.argv[1])
from partners import fincke_pohst_min

T2 = np.array([[2, 1], [1, 6]])          # trace form of O on basis (1, beta)
p = int(sys.argv[2]); mode = sys.argv[3]; ntr = int(sys.argv[4]); out = sys.argv[5]
roots = sorted(r for r in range(p) if (r * r - r + 3) % p == 0)
assert len(roots) == 2, f"{p} not split"
inv = pow((roots[1] - roots[0]) % p, -1, p)

def lift(r1, r2):
    b = ((r2 - r1) * inv) % p
    a = (r1 - b * roots[0]) % p
    return a, b

def perp_basis(basis):
    if not basis:
        return [[1 if i == j else 0 for i in range(4)] for j in range(4)]
    M = flint.nmod_mat([[int(x) for x in b] for b in basis], p)
    X, nul = M.nullspace()
    return [[int(X[i, j]) for i in range(4)] for j in range(nul)]

def build(C1b):
    C2b = perp_basis(C1b)
    gens = []
    for i in range(4):
        for e in ((p, 0), (0, p)):
            v = [0] * 8; v[2 * i], v[2 * i + 1] = e; gens.append(v)
    for c in C1b:
        v = []
        for i in range(4): v += list(lift(c[i] % p, 0))
        gens.append(v)
    for c in C2b:
        v = []
        for i in range(4): v += list(lift(0, c[i] % p))
        gens.append(v)
    Hn = flint.fmpz_mat([[int(t) for t in g] for g in gens]).hnf()
    B = np.array([[int(x) for x in r] for r in Hn.tolist()], dtype=np.int64)
    B = B[np.abs(B).sum(1) > 0]
    T8 = np.kron(np.eye(4, dtype=np.int64), T2)
    G = B @ T8 @ B.T
    assert (G % p == 0).all()
    return B, G // p, C2b

def rref_subspaces(maxdim=3):
    for k in range(1, maxdim + 1):
        for piv in itertools.combinations(range(4), k):
            free = [(r, c) for r in range(k) for c in range(piv[r] + 1, 4) if c not in piv]
            for vals in itertools.product(range(p), repeat=len(free)):
                rows = [[0] * 4 for _ in range(k)]
                for r in range(k): rows[r][piv[r]] = 1
                for (r, c), v in zip(free, vals): rows[r][c] = v
                yield rows

rng = random.Random(p)
if mode == "all":
    gen = rref_subspaces()
else:
    def rnd():
        for _ in range(ntr):
            k = rng.choice([1, 2, 2, 3])
            rows = [[rng.randrange(p) for _ in range(4)] for _ in range(k)]
            if flint.nmod_mat(rows, p).rank() == k:
                yield rows
    gen = rnd()
best = []
from collections import Counter
dist = Counter()
for C1b in gen:
    B, G, C2b = build(C1b)
    det = int(flint.fmpz_mat(G.tolist()).det())
    if det != 11 ** 4:
        continue
    mn, cnt = fincke_pohst_min(G.astype(float), 6.5)
    mn = round(mn) if mn < 1e9 else 99
    dist[mn] += 1
    if mn >= 6:
        best.append(dict(p=p, C1=C1b, C2=C2b, trace_min=mn, kissing=cnt, B=B.tolist(), G=G.tolist()))
        print(f"  FOUND trace-min {mn} kissing {cnt} C1={C1b}", flush=True)
        if len(best) >= 20:
            break
print(f"p={p} ({mode}): trace-min distribution {dict(sorted(dist.items()))}; with min>=6: {len(best)}", flush=True)
json.dump(best, open(out, "w"))
