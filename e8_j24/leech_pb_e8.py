"""Leech as Barnes P_b (x) E8 over O = Z[alpha], alpha = (1+sqrt(-7))/2.
P_b: basis (1,1,alpha), (0,beta,beta), (0,0,2), beta = 1 - alpha, h = 1/2 sum a_i conj(b_i).
L = span{ (p1 x, p2 x, p3 x) } in E8^3, Z-form = 1/2 sum b_E8(y_i, y_i)."""
import sys, json
import numpy as np
import flint
sys.path.insert(0, "/home/claude/gamma72-leech96/leech96_hermitian")
from svp_util import shortest_vector, count_short
from fpylll import IntegerMatrix
e8 = json.load(open(sys.argv[1]))[0]
F = np.array(e8["F"], dtype=np.int64); A = np.array(e8["A"], dtype=np.int64); B = np.array(e8["B"], dtype=np.int64)
I = np.eye(8, dtype=np.int64)
M = lambda a, b: a * I + b * A          # matrix of a + b*alpha
Pb = [[(1, 0), (1, 0), (0, 1)], [(0, 0), (1, -1), (1, -1)], [(0, 0), (0, 0), (2, 0)]]
gens = []
for row in Pb:
    for mult in ((1, 0), (0, 1)):        # row and alpha*row
        mats = []
        for (a, b) in row:
            # (a + b al) * mult ; al^2 = al - 2
            c, d = mult
            prod = (a * c - 2 * b * d, a * d + b * c + b * d)
            mats.append(M(*prod))
        for t in range(8):
            gens.append(np.concatenate([m[t] for m in mats]))
Hn = flint.fmpz_mat([[int(v) for v in g] for g in gens]).hnf()
Y = np.array([[int(v) for v in r] for r in Hn.tolist()], dtype=np.int64)
Y = Y[np.abs(Y).sum(1) > 0]
assert Y.shape[0] == 24
F3 = np.kron(np.eye(3, dtype=np.int64), F)
G2 = Y @ F3 @ Y.T
assert (G2 % 2 == 0).all()
G = G2 // 2
det = int(flint.fmpz_mat(G.tolist()).det()); even = (np.diag(G) % 2 == 0).all()
amb = Y @ np.kron(np.eye(3, dtype=np.int64), B)          # norm = |amb|^2 / 4
n2, v = shortest_vector(IntegerMatrix.from_matrix(amb.tolist()))
cnt = len(count_short(IntegerMatrix.from_matrix(amb.tolist()), 16, max_sols=200000))
print(f"det {det}, even {even}, minimum {n2/4}, vectors of norm 4 (up to sign): {cnt}")
# alpha action in L coordinates
A3 = np.kron(np.eye(3, dtype=np.int64), A)
AL = flint.fmpq_mat((Y @ A3).tolist()) * flint.fmpq_mat(Y.tolist()).inv()
assert all(int(e.q) == 1 for e in AL.entries())
AL = np.array([[int(str(AL[i, j])) for j in range(24)] for i in range(24)], dtype=np.int64)
I24 = np.eye(24, dtype=np.int64)
assert ((AL @ AL) - AL + 2 * I24 == 0).all() and (AL @ G @ AL.T == 2 * G).all()
json.dump(dict(G=G.tolist(), Y=Y.tolist(), AL=AL.tolist(), amb=amb.tolist()), open(sys.argv[2], "w"))
print("saved Leech with O-structure")
