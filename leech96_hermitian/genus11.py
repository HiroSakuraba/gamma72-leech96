"""Random Kneser-neighbour walk in the genus of rank-4 Hermitian unimodular lattices over
O = Z[beta], beta = (1+sqrt(-11))/2, beta^2 = beta - 3, starting from O^4 with h = sum x_i conj(y_i).
Neighbour at the principal prime (beta) of norm 3 (conj(beta) = 1 - beta):
   P_v = { x in P : h(x, v) in beta*O },   P' = P_v + O * v / conj(beta)     (v in P, 3 | h(v,v), v not in beta P)
Records the trace lattices (Tr h, 8-dim, det 11^4) and their minima; target: trace min 6 (Hermitian 3)."""
import sys, json, random, math
import numpy as np
from fractions import Fraction as Fr
import flint
sys.path.insert(0, sys.argv[1])
from partners import fincke_pohst_min

steps = int(sys.argv[2]); seed = int(sys.argv[3]); out = sys.argv[4]
rng = random.Random(seed)

# K elements: (a, b) = a + b*beta, beta^2 = beta - 3
def mul(s, t):
    a, b = s; c, d = t
    return (a * c - 3 * b * d, a * d + b * c + b * d)
def conj(s):
    a, b = s
    return (a + b, -b)
def add(s, t):
    return (s[0] + t[0], s[1] + t[1])
def herm(x, y):  # sum x_i conj(y_i)
    r = (Fr(0), Fr(0))
    for xi, yi in zip(x, y):
        r = add(r, mul(xi, conj(yi)))
    return r
BETA = (Fr(0), Fr(1)); BETAB = conj(BETA)          # 1 - beta
NBETAB_INV = (Fr(1, 3) * BETA[0], Fr(1, 3) * BETA[1])  # 1/conj(beta) = beta / 3

def to_vec8(x):
    return [c for xi in x for c in xi]
def from_vec8(v):
    return [(v[2 * i], v[2 * i + 1]) for i in range(4)]

def zbasis(gens):
    """Z-basis (8 vectors) of the Z-span of rational 8-vectors."""
    den = 1
    for g in gens:
        for c in g:
            den = math.lcm(den, Fr(c).denominator)
    M = flint.fmpz_mat([[int(Fr(c) * den) for c in g] for g in gens]).hnf()
    rows = [[Fr(int(t), den) for t in r] for r in M.tolist() if any(int(t) for t in r)]
    assert len(rows) == 8
    return rows

def trace_gram(B):
    vs = [from_vec8(b) for b in B]
    G = [[2 * herm(vs[i], vs[j])[0] + herm(vs[i], vs[j])[1] for j in range(8)] for i in range(8)]  # Tr(a+b beta) = 2a + b
    return G

def oclose(vecs):
    """Add beta-multiples so Z-span is an O-module."""
    out = []
    for v in vecs:
        x = from_vec8(v)
        out.append(v)
        out.append(to_vec8([mul(BETA, xi) for xi in x]))
    return out

def neighbour(B):
    vs = [from_vec8(b) for b in B]
    for _ in range(200):
        coeffs = [rng.randrange(3) for _ in range(8)]
        if not any(coeffs):
            continue
        v = [(Fr(0), Fr(0))] * 4
        for c, b in zip(coeffs, vs):
            v = [add(vi, (c * bi[0], c * bi[1])) for vi, bi in zip(v, b)]
        n = herm(v, v)
        if n[1] != 0 or n[0].denominator != 1 or int(n[0]) % 3 != 0:
            continue
        # f(x) = h(x, v) mod (beta): a + b beta -> a mod 3
        f = []
        for b in vs:
            s = herm(b, v)
            assert s[0].denominator == 1 and s[1].denominator == 1
            f.append(int(s[0]) % 3)
        if all(x == 0 for x in f):
            continue  # v in conj stuff, skip
        # kernel of f mod 3 on the Z-basis
        j0 = next(j for j in range(8) if f[j])
        inv = pow(f[j0], -1, 3)
        gens = []
        for j in range(8):
            if j == j0:
                gens.append([3 * c for c in B[j]])
            else:
                t = (-f[j] * inv) % 3
                gens.append([B[j][k] + t * B[j0][k] for k in range(8)])
        w = [mul(vi, NBETAB_INV) for vi in v]
        gens.append(to_vec8(w))
        gens = oclose(gens)
        try:
            NB = zbasis(gens)
        except AssertionError:
            continue
        return NB
    return None

B = [[Fr(1 if k == j else 0) for k in range(8)] for j in range(8)]  # O^4 with basis e_i, beta e_i
seen = {}
best = []
for st in range(steps):
    NB = neighbour(B)
    if NB is None:
        continue
    G = trace_gram(NB)
    if any(Fr(g).denominator != 1 for r in G for g in r):
        continue
    from gram_lll import gram_lll
    Uarr, Gred = gram_lll([[int(g) for g in r] for r in G])
    NB = [[sum(Uarr[i][k] * NB[k][c] for k in range(8)) for c in range(8)] for i in range(8)]
    Gi = np.array([[int(g) for g in r] for r in trace_gram(NB)], dtype=np.int64)
    assert (Gi == np.array(Gred)).all()
    det = int(flint.fmpz_mat(Gi.tolist()).det())
    if det != 11 ** 4:
        continue
    mn, cnt = fincke_pohst_min(Gi.astype(float), 6.5)
    mn = round(mn) if mn < 1e9 else 99
    key = (mn, cnt)
    seen[key] = seen.get(key, 0) + 1
    B = NB  # walk
    if mn >= 6:
        best.append(dict(trace_min=mn, kissing=cnt, B=[[str(c) for c in r] for r in NB], G=Gi.tolist()))
        print(f"step {st}: FOUND trace-min {mn}, kissing {cnt}", flush=True)
        if len(best) >= 5:
            break
    if st % 200 == 0:
        print(f"step {st}: (min, kissing) classes seen {sorted(seen.items())}", flush=True)
print("final classes seen (trace-min, kissing): count", sorted(seen.items()))
json.dump(best, open(out, "w"))
