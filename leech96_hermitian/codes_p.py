"""Search length-4 codes C1 over F_p (C2 = C1^perp) whose Euclidean (Lee-squared) minimum weights
clear the twist thresholds, build the partner P_C, check Hermitian unimodularity and min(P) >= 3.

Euclidean weight of c in F_p^n: sum_i min(c_i, p - c_i)^2.
Threshold: words of C1 meet the twist by prime #2 (root r2), words of C2 the twist by prime #1.
need  ew(C1) * m2 >= 10  and  ew(C2) * m1 >= 10, m = rescaled min of the twist (2 or 4)."""
import sys, json, itertools, random, math
import numpy as np
import flint
sys.path.insert(0, sys.argv[1])
from partners import fincke_pohst_min, T2

D = sys.argv[1]
ideal_min = json.load(open(f"{D}/ideal_minima.json"))
sign = int(sys.argv[2]); p = int(sys.argv[3]); trials = int(sys.argv[4]); out = sys.argv[5]
rng = random.Random(p * 7 + sign)
roots = sorted(r for r in range(p) if (r * r - r + 6) % p == 0)
m1 = ideal_min[f"sign{sign}_p{p}_r{roots[0]}"]["min_over_p"]
m2 = ideal_min[f"sign{sign}_p{p}_r{roots[1]}"]["min_over_p"]
need1 = math.ceil(10 / m2 - 1e-9)   # min Euclidean weight for C1
need2 = math.ceil(10 / m1 - 1e-9)   # for C2
print(f"sign {sign} p={p}: twist minima m1={m1} m2={m2}; need ew(C1)>={need1}, ew(C2)>={need2}", flush=True)

def ew(v):
    return sum(min(x % p, p - x % p) ** 2 for x in v)

def span(basis):
    for coeffs in itertools.product(range(p), repeat=len(basis)):
        if any(coeffs):
            yield [sum(c * b[i] for c, b in zip(coeffs, basis)) % p for i in range(4)]

def min_ew(basis):
    return min(ew(v) for v in span(basis)) if basis else math.inf

def perp_basis(basis):
    k = len(basis)
    M = flint.nmod_mat([[int(x) for x in b] for b in basis], p) if k else None
    if k == 0:
        return [[1 if i == j else 0 for i in range(4)] for j in range(4)]
    ns = M.nullspace()[0]  # columns span kernel
    rank = 4 - M.rank()
    cols = [[int(ns[i, j]) for i in range(4)] for j in range(rank)]
    return cols

def lift(r1, r2):
    inv = pow((roots[1] - roots[0]) % p, -1, p)
    b = ((r2 - r1) * inv) % p
    a = (r1 - b * roots[0]) % p
    return a, b

def build_P(C1b, C2b):
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
    Hn = flint.fmpz_mat(gens).hnf()
    B = np.array([[int(x) for x in r] for r in Hn.tolist()], dtype=np.int64)
    B = B[np.abs(B).sum(1) > 0]
    T8 = np.kron(np.eye(4, dtype=np.int64), T2)
    G = B @ T8 @ B.T
    assert (G % p == 0).all()
    G //= p
    return B, G

found = []
seen = set()
for t in range(trials):
    k = rng.choice([1, 2, 3])
    C1b = [[rng.randrange(p) for _ in range(4)] for _ in range(k)]
    if flint.nmod_mat(C1b, p).rank() < k:
        continue
    key = tuple(sorted(tuple(v) for v in span(C1b))) if p ** k < 3000 else None
    e1 = min_ew(C1b)
    if e1 < need1:
        continue
    C2b = perp_basis(C1b)
    e2 = min_ew(C2b)
    if e2 < need2:
        continue
    B, G = build_P(C1b, C2b)
    det = int(flint.fmpz_mat(G.tolist()).det())
    if det != 23 ** 4:
        continue
    mn, cnt = fincke_pohst_min(G.astype(float), 8.5)
    if mn < 6 - 1e-6:
        continue
    found.append(dict(p=p, sign=sign, C1=C1b, C2=C2b, ew1=e1, ew2=e2, trace_min=round(mn), kissing=cnt,
                      B=B.tolist(), G=G.tolist()))
    print(f"  k={k} ew(C1)={e1} ew(C2)={e2} trace-min(P)={round(mn)} kiss={cnt}", flush=True)
    if len(found) >= 6:
        break
json.dump(found, open(out, "w"))
print("found", len(found))
