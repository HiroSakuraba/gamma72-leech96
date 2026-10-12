"""Symmetric search for an extremal 7-modular J inside Leech = P_b (x) E8 over Z[(1+sqrt-7)/2].
Automorphisms: mu (x -> 2x on the F_7 coordinates of E8, order 3, commutes with sqrt(-7)) applied to all
three E8 blocks, optionally combined with signed coordinate permutations of P_b that preserve P_b.
For g of order 3 acting on V = Leech/sqrt(-7)Leech (F_7^12): V = V_1 + V_2 + V_4 (eigenvalues = cube roots
of 1 in F_7). psi pairs V_lambda with V_(1/lambda): V_1 with itself, V_2 with V_4.  g-invariant Lagrangians:
C = L_1 + W + (W^perp in V_4), L_1 Lagrangian in V_1, W any subspace of V_2.  Sample them; score = number of
vector pairs of norm <= 8 in M_C (0 = extremal)."""
import sys, json, random, itertools, time
import numpy as np
import flint
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "leech96_hermitian"))
from svp_util import count_short
from fpylll import IntegerMatrix
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_decomp
e8 = json.load(open("e8_seven.json"))[0]; D = json.load(open("leech_o7.json"))
samples = int(sys.argv[1]); seed = int(sys.argv[2]); out = sys.argv[3]
rng = random.Random(seed)
F8 = np.array(e8["F"]); A8 = np.array(e8["A"]); B8 = np.array(e8["B"])
G = np.array(D["G"], dtype=np.int64); AL = np.array(D["AL"], dtype=np.int64); amb = np.array(D["amb"], dtype=np.int64)
Y = np.array(D["Y"], dtype=np.int64)
I24 = np.eye(24, dtype=np.int64); S = 2 * AL - I24
# --- E8 automorphism mu: coordinate i -> 2i on F_7, infinity fixed
P = np.zeros((8, 8), dtype=np.int64)
for i in range(7): P[i, (2 * i) % 7] = 1
P[7, 7] = 1
Bq = flint.fmpq_mat(B8.tolist())
X = Bq * flint.fmpq_mat(P.tolist()) * Bq.inv()
assert all(int(e.q) == 1 for e in X.entries())
mu = np.array([[int(str(X[i, j])) for j in range(8)] for i in range(8)], dtype=np.int64)
assert (mu @ A8 == A8 @ mu).all() and (mu @ F8 @ mu.T == F8).all()
# --- P_b automorphisms among signed permutations of O^3 (rows as O-vectors), O = Z[al], al^2 = al - 2
Pb = [[(1, 0), (1, 0), (0, 1)], [(0, 0), (1, -1), (1, -1)], [(0, 0), (0, 0), (2, 0)]]
def M(a, b): return a * np.eye(8, dtype=np.int64) + b * A8
def pb_gens(perm, signs):
    gens = []
    for row in Pb:
        r2 = [None] * 3
        for k in range(3):
            a, b = row[k]; r2[perm[k]] = (signs[k] * a, signs[k] * b)
        for mult in ((1, 0), (0, 1)):
            mats = []
            for (a, b) in r2:
                c, d = mult; mats.append(M(a * c - 2 * b * d, a * d + b * c + b * d))
            for t in range(8):
                gens.append(np.concatenate([m[t] for m in mats]))
    return gens
def same_lattice(gens):
    Hn = flint.fmpz_mat([[int(v) for v in g] for g in gens]).hnf()
    Z = np.array([[int(v) for v in r] for r in Hn.tolist()], dtype=np.int64)
    Z = Z[np.abs(Z).sum(1) > 0]
    return Z.shape == Y.shape and (Z == Y).all()
pb_auts = []
for perm in itertools.permutations(range(3)):
    for signs in itertools.product([1, -1], repeat=3):
        if same_lattice(pb_gens(perm, signs)):
            pb_auts.append((perm, signs))
print("signed permutations preserving P_b:", len(pb_auts), flush=True)
def leech_aut(perm, signs, use_mu):
    """matrix on E8^3 coordinates (row vectors): block k goes to block perm[k] with sign, then mu."""
    Mbig = np.zeros((24, 24), dtype=np.int64)
    m = mu if use_mu else np.eye(8, dtype=np.int64)
    for k in range(3):
        Mbig[8 * k:8 * k + 8, 8 * perm[k]:8 * perm[k] + 8] = signs[k] * m
    Xl = flint.fmpq_mat((Y @ Mbig).tolist()) * flint.fmpq_mat(Y.tolist()).inv()
    if not all(int(e.q) == 1 for e in Xl.entries()): return None
    gL = np.array([[int(str(Xl[i, j])) for j in range(24)] for i in range(24)], dtype=np.int64)
    assert (gL @ G @ gL.T == G).all() and (gL @ AL == AL @ gL).all()
    return gL
# --- V and psi
Dg, U, V = smith_normal_decomp(Matrix(S.tolist()))
V = np.array(V.tolist(), dtype=object); Vinv = np.array(Matrix(V.tolist()).inv().tolist(), dtype=object)
idx7 = [i for i in range(24) if abs(int(Dg[i, i])) == 7]
def proj(c):
    w = np.array(c, dtype=object).dot(V); return np.array([int(w[i]) % 7 for i in idx7], dtype=np.int64)
lifts = []
for i in idx7:
    e = np.zeros(24, dtype=object); e[i] = 1; lifts.append(np.array([int(x) for x in e.dot(Vinv)], dtype=object))
SG = S.astype(object).dot(G.astype(object))
Psi = np.array([[int(lifts[i].dot(SG).dot(lifts[j])) % 7 for j in range(12)] for i in range(12)], dtype=np.int64)
def vmat(gL):
    return np.array([proj(lifts[i].dot(gL.astype(object))) for i in range(12)], dtype=np.int64)
def nullsp(Mx):
    X, k = flint.nmod_mat([[int(v) for v in r] for r in Mx], 7).nullspace()
    return [np.array([int(X[i, j]) for i in range(Mx.shape[1])], dtype=np.int64) for j in range(k)]
def rowspace_basis(rows):
    if not rows: return []
    Mx = flint.nmod_mat([[int(v) for v in r] for r in rows], 7).rref()[0]
    return [np.array([int(Mx[i, j]) for j in range(12)]) for i in range(Mx.nrows()) if any(int(Mx[i, j]) for j in range(12))]
def eigenspace(gv, lam):
    # row vectors x with x gv = lam x  <=>  (gv^T - lam I) x^T = 0
    return nullsp(((gv.T - lam * np.eye(12, dtype=np.int64)) % 7))
def lattice_of(C):
    gens = [list(map(int, r)) for r in S] + [list(map(int, sum(int(C[k][i]) * lifts[i] for i in range(12)))) for k in range(len(C))]
    Hn = flint.fmpz_mat(gens).hnf()
    Z = np.array([[int(v) for v in r] for r in Hn.tolist()], dtype=np.int64)
    return Z[np.abs(Z).sum(1) > 0]
def score(C):
    Z = lattice_of(C); E = Z @ amb
    vs = count_short(IntegerMatrix.from_matrix(E.tolist()), 32, max_sols=100000)
    ns = sorted(sum(x * x for x in v) // 4 for v in vs)
    return len(vs), (ns[0] if ns else None), Z
def rand_vec_in(basis):
    return sum(rng.randrange(7) * b for b in basis) % 7 if basis else None
def random_lagrangian_in(basis):
    """random Lagrangian subspace of the psi-nondegenerate subspace spanned by basis"""
    out = []
    target = len(basis) // 2
    tries = 0
    while len(out) < target and tries < 5000:
        tries += 1
        v = rand_vec_in(basis)
        if not v.any(): continue
        if all((v @ Psi @ b) % 7 == 0 for b in out) and len(rowspace_basis(out + [v])) == len(out) + 1:
            out.append(v)
    return out if len(out) == target else None
groups = []
for (perm, signs) in pb_auts:
    for use_mu in (True,):
        gL = leech_aut(perm, signs, use_mu)
        if gL is None: continue
        gv = vmat(gL)
        g3 = (gv @ gv @ gv) % 7
        if not (g3 == np.eye(12, dtype=np.int64)).all(): continue   # need order dividing 3 on V
        E1, E2, E4 = eigenspace(gv, 1), eigenspace(gv, 2), eigenspace(gv, 4)
        if len(E1) + len(E2) + len(E4) != 12: continue
        groups.append((perm, signs, gv, E1, E2, E4))
print("order-3 symmetries usable:", [(g[0], g[1], len(g[3]), len(g[4]), len(g[5])) for g in groups], flush=True)
best = None; t0 = time.time(); hist = {}
for s in range(samples):
    perm, signs, gv, E1, E2, E4 = groups[s % len(groups)]
    L1 = random_lagrangian_in(E1) if E1 else []
    if L1 is None: continue
    k = rng.randrange(0, len(E2) + 1)
    W = rowspace_basis([rand_vec_in(E2) for _ in range(k)]) if k else []
    # W^perp inside E4: x in E4 with x Psi w = 0 for all w in W
    if W:
        coeffs = nullsp(np.array([[int((e @ Psi @ w) % 7) for e in E4] for w in W]).T.T) if False else None
        Mx = np.array([[int((e @ Psi @ w) % 7) for e in E4] for w in W])
        ker = nullsp(Mx)
        Wp = [sum(int(c[i]) * E4[i] for i in range(len(E4))) % 7 for c in ker]
    else:
        Wp = list(E4)
    C = rowspace_basis(list(L1) + list(W) + list(Wp))
    if len(C) != 6: continue
    assert all((a @ Psi @ b) % 7 == 0 for a in C for b in C)
    sc, mn, Z = score(C)
    hist[sc] = hist.get(sc, 0) + 1
    if best is None or sc < best[0]:
        best = (sc, mn, [list(map(int, c)) for c in C], perm, signs)
        print(f"sample {s}: {sc} short pairs (min {mn}) symmetry perm={perm} signs={signs} [{time.time()-t0:.0f}s]", flush=True)
        json.dump(dict(short_pairs=sc, min=mn, C=best[2], perm=list(perm), signs=list(signs)), open(out, "w"))
        if sc == 0:
            print("FOUND extremal J", flush=True); break
print("best", best[:2], "distribution of scores (lowest 10):", sorted(hist.items())[:10])
