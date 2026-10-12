"""Build L = P (x)_O Leech  (96-dim) for a partner P from partners*.json and a Leech O-structure,
verify even unimodular exactly, then try to find vectors of norm <= 8 with LLL/BKZ (kill test).

Ambient integer embedding: y in Z^96 (four Leech ambient blocks, Leech scaled by sqrt 8),
norm(y) = |y|^2 / (8 p) with p the partner prime (p = 3 here)."""
import sys, json, time
import numpy as np
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
from fpylll import IntegerMatrix, LLL, BKZ, GSO
from fpylll.algorithms.bkz2 import BKZReduction


def build(P, st):
    F = np.array(st["F"], dtype=np.int64)
    A = np.array(st["A"], dtype=np.int64)
    H = np.array(st["H"], dtype=np.int64)
    I = np.eye(24, dtype=np.int64)
    # sanity: alpha^2 - alpha + 6 = 0 ; A F A^T = 6F ; adjoint relation A F + F A^T = F
    assert ((A @ A) - A + 6 * I == 0).all()
    assert (A @ F @ A.T == 6 * F).all()
    assert (A @ F + F @ A.T == F).all()
    p = P["p"]
    Bp = np.array(P["B"], dtype=np.int64)          # 8 rows: (a1,b1,...,a4,b4)
    gens = []
    for row in Bp:
        Ms = [row[2 * i] * I + row[2 * i + 1] * A for i in range(4)]
        for t in range(24):
            gens.append(np.concatenate([Ms[i][t] for i in range(4)]))
    gens = np.array(gens)
    import flint
    Hn = flint.fmpz_mat([[int(v) for v in r] for r in gens]).hnf()
    Y = np.array([[int(v) for v in r] for r in Hn.tolist()], dtype=object)
    Y = Y[[any(v != 0 for v in r) for r in Y]]
    assert Y.shape[0] == 96, Y.shape
    F4 = np.kron(np.eye(4, dtype=np.int64), F).astype(object)
    Gm = Y.dot(F4).dot(Y.T)
    assert all(int(v) % p == 0 for v in Gm.flatten()), "not integral before scaling"
    Gm = Gm // p
    det = int(flint.fmpz_mat([[int(v) for v in r] for r in Gm]).det())
    even = all(int(Gm[i, i]) % 2 == 0 for i in range(96))
    H4 = np.kron(np.eye(4, dtype=np.int64), H).astype(object)
    Yamb = Y.dot(H4)
    return Y, Gm, det, even, Yamb, p

