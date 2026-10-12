"""Search for an extremal 7-modular lattice J (dim 24, min 10) as M_C = pi^{-1}(C) inside Leech,
C a Lagrangian subspace of V = Leech / sqrt(-7) Leech (12-dim over F_7) for the symplectic form
psi(z, m) = b(z S, m) mod 7, S = 2 alpha - 1.  M_C^# = (1/7) S M_C, so M_C is 7-modular, det 7^12.
Hill-climb with random symplectic transvections; score = number of vector pairs of norm <= 8."""
import sys, json, random, time
import numpy as np
import flint
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "leech96_hermitian"))
from svp_util import count_short
from fpylll import IntegerMatrix
D = json.load(open(sys.argv[1])); iters = int(sys.argv[2]); seed = int(sys.argv[3]); out = sys.argv[4]
rng = random.Random(seed)
G = np.array(D["G"], dtype=np.int64); AL = np.array(D["AL"], dtype=np.int64); amb = np.array(D["amb"], dtype=np.int64)
I = np.eye(24, dtype=np.int64); S = 2 * AL - I
assert ((S @ S) == -7 * I).all()
# Smith form of S: U S V = diag  (flint returns snf only; compute via sympy for transforms)
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_decomp
Dg, U, V = smith_normal_decomp(Matrix(S.tolist()))
Dg = np.array(Dg.tolist(), dtype=object); U = np.array(U.tolist(), dtype=object); V = np.array(V.tolist(), dtype=object)
idx7 = [i for i in range(24) if abs(int(Dg[i, i])) == 7]
assert len(idx7) == 12 and all(abs(int(Dg[i, i])) in (1, 7) for i in range(24))
Vinv = np.array(Matrix(V.tolist()).inv().tolist(), dtype=object)
# row-vector convention: sublattice S L = rows of S.  c in SL  <=>  (c V)_i divisible by d_i  (since S = U^-1 D V^-1)
def proj(c):
    w = np.array(c, dtype=object).dot(V)
    return np.array([int(w[i]) % 7 for i in idx7], dtype=np.int64)
lifts = []
for i in idx7:
    e = np.zeros(24, dtype=object); e[i] = 1
    lifts.append(np.array([int(x) for x in e.dot(Vinv)], dtype=object))
for k, l in enumerate(lifts):
    pv = proj(l); assert pv[k] == 1 and pv.sum() % 7 == 1
SG = S.astype(object).dot(G.astype(object))
Psi = np.array([[int(lifts[i].dot(SG).dot(lifts[j])) % 7 for j in range(12)] for i in range(12)], dtype=np.int64)
assert ((Psi + Psi.T) % 7 == 0).all() and int(flint.nmod_mat(Psi.tolist(), 7).rank()) == 12
def lagrangian_random():
    # greedy random isotropic basis
    basis = []
    while len(basis) < 6:
        v = np.array([rng.randrange(7) for _ in range(12)], dtype=np.int64)
        if not v.any(): continue
        if all((v @ Psi @ b) % 7 == 0 for b in basis) and \
           flint.nmod_mat([list(map(int, b)) for b in basis] + [list(map(int, v))], 7).rank() == len(basis) + 1:
            basis.append(v)
    return np.array(basis)
def transvect(C):
    v = np.array([rng.randrange(7) for _ in range(12)], dtype=np.int64)
    c = rng.randrange(1, 7)
    # tau(x) = x + c psi(x, v) v
    return np.array([(x + c * int(x @ Psi @ v) * v) % 7 for x in C])
def lattice_of(C):
    gens = [list(map(int, r)) for r in S] + [list(map(int, sum(int(C[k, i]) * lifts[i] for i in range(12)))) for k in range(6)]
    Hn = flint.fmpz_mat(gens).hnf()
    Y = np.array([[int(v) for v in r] for r in Hn.tolist()], dtype=np.int64)
    return Y[np.abs(Y).sum(1) > 0]
def score(C, check=False):
    Y = lattice_of(C)
    assert Y.shape[0] == 24
    if check:
        GM = Y @ G @ Y.T
        assert int(flint.fmpz_mat(GM.tolist()).det()) == 7 ** 12
    E = Y @ amb
    vs = count_short(IntegerMatrix.from_matrix(E.tolist()), 32, max_sols=100000)
    norms = sorted(sum(x * x for x in v) // 4 for v in vs)
    return len(vs), (norms[0] if norms else None)
C = lagrangian_random()
cur, mn = score(C, check=True)
best = (cur, mn, C.copy())
t0 = time.time()
print(f"start: {cur} short pairs, min {mn}", flush=True)
for it in range(iters):
    C2 = transvect(C)
    if flint.nmod_mat([list(map(int, r)) for r in C2], 7).rank() < 6: continue
    s2, m2 = score(C2)
    if s2 <= cur:
        C, cur, mn = C2, s2, m2
        if cur < best[0]:
            best = (cur, mn, C.copy())
            print(f"it {it}: {cur} short pairs (min {mn})  [{time.time()-t0:.0f}s]", flush=True)
            json.dump(dict(short_pairs=cur, min=mn, C=C.tolist()), open(out, "w"))
            if cur == 0:
                print("FOUND extremal J (no vectors of norm <= 8)", flush=True); break
print("best", best[0], "min", best[1])
