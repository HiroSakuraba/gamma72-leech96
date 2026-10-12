"""Z[beta]-structures on Leech, beta = (1+sqrt(-11))/2 (beta^2 - beta + 3 = 0), from the order-11
automorphism tau: i -> 4i on F_23 coordinates (fixes 0 and infinity) of the QR-Golay Leech lattice.
sqrt(-11) = g + T where g is the Gauss sum over <tau> (acts on the 20-dim non-fixed part) and
T is a sqrt(-11) on the 4-dim fixed sublattice, found by exhaustive search."""
import sys, json, itertools, math
import numpy as np
import flint

D = sys.argv[1]
st23 = json.load(open(f"{D}/leech23.json"))[0]
F = np.array(st23["F"], dtype=np.int64); H = np.array(st23["H"], dtype=np.int64)   # Leech basis rows (ambient, sqrt8-scaled)
Hq = flint.fmpq_mat(H.tolist()); Hinv = Hq.inv()
p = 23
# tau permutation (row-vector action x -> x Ptau): coordinate i -> 4i mod 23, infinity (index 23) fixed
Pt = np.zeros((24, 24), dtype=np.int64)
for i in range(p):
    Pt[i, (4 * i) % p] = 1
Pt[23, 23] = 1

def leech_coords(Mamb):
    X = Hq * flint.fmpq_mat([[int(v) for v in r] for r in Mamb]) * Hinv
    if all(int(e.q) == 1 for e in X.entries()):
        return np.array([[int(str(X[i, j])) for j in range(24)] for i in range(24)], dtype=np.int64)
    return None

assert leech_coords(Pt) is not None, "tau not an automorphism"
leg11 = {0: 0}
for k in range(1, 11):
    leg11[k] = 1 if pow(k, 5, 11) == 1 else -1
g = np.zeros((24, 24), dtype=np.int64)
Pk = np.eye(24, dtype=np.int64)
for k in range(11):
    if leg11[k]:
        g += leg11[k] * Pk
    Pk = Pk @ Pt
# fixed sublattice of tau in Leech (saturated): LLL on [W*(tauL - I) | I]
from fpylll import IntegerMatrix, LLL
tauL = leech_coords(Pt)
Mk = tauL - np.eye(24, dtype=np.int64)
W = 10**6
aug = np.hstack([W * Mk, np.eye(24, dtype=np.int64)])
AM = IntegerMatrix.from_matrix(aug.tolist())
LLL.reduction(AM)
rows = [[AM[i, j] for j in range(48)] for i in range(24)]
Bf = np.array([r[24:] for r in rows if all(x == 0 for x in r[:24])], dtype=np.int64)
assert (Bf @ tauL == Bf).all()
snf = flint.fmpz_mat(Bf.tolist()).snf()
assert all(abs(int(snf[i, i])) == 1 for i in range(Bf.shape[0])), "not saturated"
Ured = IntegerMatrix.identity(Bf.shape[0])
Mred = IntegerMatrix.from_matrix((Bf @ H).tolist())
LLL.reduction(Mred, Ured)
Uarr = np.array([[Ured[i, j] for j in range(Bf.shape[0])] for i in range(Bf.shape[0])], dtype=np.int64)
Bf = Uarr @ Bf
Gf = Bf @ F @ Bf.T
print("fixed sublattice rank", Bf.shape[0], "Gram", Gf.tolist(), "det", int(round(np.linalg.det(Gf))), flush=True)

# search T (4x4 integer, rows in fixed-lattice coords) with T Gf T^T = 11 Gf and T^2 = -11
k = Bf.shape[0]
def vecs_of_norm(nrm):
    out = []
    R = 4
    for c in itertools.product(range(-R, R + 1), repeat=k):
        c = np.array(c)
        if c @ Gf @ c == nrm:
            out.append(c)
    return out
cands = [vecs_of_norm(11 * Gf[i, i]) for i in range(k)]
print("candidate row counts", [len(c) for c in cands], flush=True)
sols = []
def rec(rows):
    i = len(rows)
    if i == k:
        T = np.array(rows)
        if ((T @ T) == -11 * np.eye(k, dtype=np.int64)).all():
            sols.append(T)
        return
    for v in cands[i]:
        if all(v @ Gf @ rows[j] == 11 * Gf[i, j] for j in range(i)):
            rec(rows + [v])
rec([])
print("T solutions on fixed lattice:", len(sols), flush=True)

# ambient: fixed lattice basis vectors in ambient coords
Bamb = Bf @ H
Bq = flint.fmpq_mat(Bamb.tolist())
Gram_amb = Bq * Bq.transpose()
Proj = Bq.transpose() * Gram_amb.inv()       # 24 x k  : x -> coefficients
results = []
for T in sols:
    Tamb = Proj * flint.fmpq_mat(T.tolist()) * Bq   # 24 x 24 rational
    S = flint.fmpq_mat(g.tolist()) + Tamb
    if not all(int(e.q) == 1 for e in S.entries()):
        # S itself may be non-integral in ambient coords yet preserve Leech; test in Leech coords
        pass
    SL = Hq * S * Hinv
    if not all(int(e.q) == 1 for e in SL.entries()):
        continue
    SLn = np.array([[int(str(SL[i, j])) for j in range(24)] for i in range(24)], dtype=np.int64)
    I = np.eye(24, dtype=np.int64)
    if not ((SLn @ SLn) == -11 * I).all():
        continue
    if not (SLn @ F @ SLn.T == 11 * F).all():
        continue
    if not ((I + SLn) % 2 == 0).all():
        continue
    A = (I + SLn) // 2
    assert ((A @ A) - A + 3 * I == 0).all()
    results.append(dict(A=A.tolist(), F=F.tolist(), H=H.tolist(), T=T.tolist()))
print("Z[beta]-structures on Leech found:", len(results), flush=True)
json.dump(results, open(f"{D}/leech11.json", "w"))
