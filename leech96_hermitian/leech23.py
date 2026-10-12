"""Leech lattice from the extended quadratic-residue code of length 24 (coords: index 0..22 = F_23,
index 23 = infinity), and a sqrt(-23) structure built from the Gauss sum over the coordinate shift.

Lattice convention: vectors x in Z^24 (scaled by sqrt 8); inner product (x.y)/8.
"""
import itertools, json, sys
import numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form

p = 23
QR = sorted({(x * x) % p for x in range(1, p)})
NQR = [x for x in range(1, p) if x not in QR]
chi = {0: 0, **{x: 1 for x in QR}, **{x: -1 for x in NQR}}


def gf2_rank_basis(rows):
    rows = [int("".join(map(str, r)), 2) for r in rows]
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r:
            basis.append(r)
    return basis


def golay():
    # try generating vectors: indicator of QR, NQR, with/without 0
    for gen in (QR, NQR, [0] + QR, [0] + NQR):
        v = [1 if i in gen else 0 for i in range(p)]
        shifts = [v[-s:] + v[:-s] for s in range(p)]
        b = gf2_rank_basis(shifts)
        if len(b) != 12:
            continue
        words = []
        for coeffs in itertools.product([0, 1], repeat=12):
            w = 0
            for c, bb in zip(coeffs, b):
                if c:
                    w ^= bb
            bits = [int(ch) for ch in format(w, f"0{p}b")]
            bits.append(sum(bits) % 2)  # parity at infinity
            words.append(bits)
        wts = sorted(set(sum(w) for w in words))
        if wts == [0, 8, 12, 16, 24]:
            return np.array(words, dtype=np.int64), gen
    raise RuntimeError("no Golay code found")


C, gen = golay()
print("Golay code from shifts of", "QR" if gen == QR else gen[:3], "weights ok, size", len(C))

# Leech generators (scaled by sqrt8): 2*c for c in C ; 4(e_i - e_j); 8 e_i ; (-3, 1^23)
gens = [2 * c for c in C[1:]]
for i in range(1, 24):
    v = np.zeros(24, dtype=np.int64); v[0] = 4; v[i] = -4; gens.append(v)
v = np.zeros(24, dtype=np.int64); v[0] = 8; gens.append(v)
v = np.ones(24, dtype=np.int64); v[0] = -3; gens.append(v)
H = hermite_normal_form(Matrix(np.array(gens).T.tolist())).T  # rows span lattice
H = np.array(H.tolist(), dtype=np.int64)
H = H[np.abs(H).sum(1) > 0]
assert H.shape[0] == 24, H.shape
F = (H @ H.T)
assert (F % 8 == 0).all()
F //= 8
from sympy import Matrix as M
print("Leech Gram det", M(F.tolist()).det(), "diag min", F.diagonal().min(), "even", (F.diagonal() % 2 == 0).all())

Hinv = np.array(M(H.tolist()).inv().tolist(), dtype=object)

def in_lattice(vec):
    c = np.array(M([list(map(int, vec))]) * M(H.tolist()).inv())
    return all(x.q == 1 for x in c.flatten())

def coords_of_ambient_map(Smat):
    """Given an ambient 24x24 integer matrix acting on row vectors (x -> x S), return its matrix
    in the Leech basis if it preserves the lattice, else None."""
    img = M(H.tolist()) * M(Smat.tolist())
    X = img * M(H.tolist()).inv()
    if all(x.q == 1 for x in X):
        return np.array(X.tolist(), dtype=np.int64)
    return None

# shift permutation sigma: coordinate i -> i+1 on F_23, fixes infinity (index 23)
P = np.zeros((24, 24), dtype=np.int64)
for i in range(p):
    P[i, (i + 1) % p] = 1
P[23, 23] = 1
assert coords_of_ambient_map(P) is not None, "shift not a Leech automorphism"
g = np.zeros((24, 24), dtype=np.int64)
Pk = np.eye(24, dtype=np.int64)
for k in range(p):
    if chi[k]:
        g += chi[k] * Pk
    Pk = Pk @ P
ones = np.zeros(24, dtype=np.int64); ones[:p] = 1
results = []
for sgn in (1, -1):
    T = np.zeros((24, 24), dtype=np.int64)
    # row-vector action: x -> x T ; T maps ones_23 -> sgn*23 e_inf, e_inf -> -sgn*ones_23
    # write x = a*ones + b*e_inf + (rest orthogonal); rest killed. Implement as rank-2 matrix:
    T += sgn * np.outer(ones, np.eye(24, dtype=np.int64)[23])          # ones-coordinate sum -> e_inf ... see check below
    T += -sgn * np.outer(np.eye(24, dtype=np.int64)[23], ones)
    S = g + T
    ok_sq = ((S @ S) == -23 * np.eye(24)).all()
    ok_sim = ((S @ S.T) == 23 * np.eye(24)).all()
    SL = coords_of_ambient_map(S)
    A_amb = None
    AL = None
    if SL is not None:
        AL2 = np.eye(24, dtype=np.int64) + SL
        if (AL2 % 2 == 0).all():
            AL = AL2 // 2
    print(f"sign {sgn}: S^2=-23 {ok_sq}, S S^T=23 {ok_sim}, preserves Leech {SL is not None}, "
          f"(1+S)/2 integral on Leech {AL is not None}")
    if AL is not None and ok_sq and ok_sim:
        results.append(dict(sign=sgn, F=F.tolist(), A=AL.tolist(), S=SL.tolist(), H=H.tolist()))
json.dump(results, open(sys.argv[1], "w"))
print("structures found:", len(results))
