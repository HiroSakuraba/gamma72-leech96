"""Kneser 2-neighbour hill-climb in the genus of 7-modular Z[alpha]-lattices (dim 24, det 7^12),
starting from the best M_C inside Leech. Lattices stored as integer coordinate rows N with a common
denominator d (row/d = vector in Leech coordinates). Step: v in M with b(v,v) = 0 mod 8, b(v, al v) = 0 mod 4,
M_v = {x : b(x,v) = b(x, al v) = 0 mod 2}, M' = M_v + Z v/2 + Z al v/2. Accepted only if M' is even,
det 7^12, and alpha-stable (exact checks). Score: vector pairs of norm <= 8 (0 = extremal)."""
import sys, json, random, time
import numpy as np
import flint
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "leech96_hermitian"))
from svp_util import count_short
from fpylll import IntegerMatrix
D = json.load(open("leech_o7.json")); start = json.load(open(sys.argv[1]))
steps = int(sys.argv[2]); seed = int(sys.argv[3]); out = sys.argv[4]
rng = random.Random(seed)
G = np.array(D["G"], dtype=object); AL = np.array(D["AL"], dtype=object); amb = np.array(D["amb"], dtype=object)
exec(open("j24_symmetric.py").read().split("# --- V and psi")[1].split("groups = []")[0]
     .replace("rng = random.Random(seed)", ""), globals()) if False else None
# rebuild starting lattice from its code C
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_decomp
S = 2 * AL - np.eye(24, dtype=object)
Dg, U, V = smith_normal_decomp(Matrix(S.tolist()))
Vinv = np.array(Matrix(V.tolist()).inv().tolist(), dtype=object)
idx7 = [i for i in range(24) if abs(int(Dg[i, i])) == 7]
lifts = []
for i in idx7:
    e = np.zeros(24, dtype=object); e[i] = 1; lifts.append(e.dot(Vinv))
C = start["C"]
gens = [list(map(int, r)) for r in S] + [list(map(int, sum(int(c[i]) * lifts[i] for i in range(12)))) for c in C]
Hn = flint.fmpz_mat(gens).hnf()
N = np.array([[int(v) for v in r] for r in Hn.tolist() if any(int(v) for v in r)], dtype=object); d = 1
def gram(N, d):
    return N.dot(G).dot(N.T)            # = d^2 * true Gram
def valid(N, d):
    Gm = gram(N, d)
    if any(int(x) % (d * d) for x in Gm.flatten()): return False
    Gt = Gm // (d * d)
    if any(int(Gt[i, i]) % 2 for i in range(24)): return False
    if int(flint.fmpz_mat([[int(x) for x in r] for r in Gt]).det()) != 7 ** 12: return False
    # alpha-stability: rows of N AL lie in Z-span of N
    X = flint.fmpq_mat([[int(x) for x in r] for r in N.dot(AL)]) * flint.fmpq_mat([[int(x) for x in r] for r in N]).inv()
    return all(int(e.q) == 1 for e in X.entries())
def score(N, d):
    E = N.dot(amb)                       # norm = |E|^2 / (4 d^2)
    vs = count_short(IntegerMatrix.from_matrix([[int(x) for x in r] for r in E]), 32 * d * d, max_sols=100000)
    ns = sorted(sum(int(x) ** 2 for x in v) // (4 * d * d) for v in vs)
    return len(vs), (ns[0] if ns else None)
def reduce_basis(N):
    E = IntegerMatrix.from_matrix([[int(x) for x in r] for r in N.dot(amb)])
    Ui = IntegerMatrix.identity(24)
    from fpylll import LLL
    LLL.reduction(E, Ui)
    Um = np.array([[Ui[i, j] for j in range(24)] for i in range(24)], dtype=object)
    return Um.dot(N)
assert valid(N, d)
cur, mn = score(N, d)
print(f"start: {cur} short pairs (min {mn})", flush=True)
best = cur; t0 = time.time(); acc = 0
for st in range(steps):
    c = np.array([rng.randrange(-1, 2) for _ in range(24)], dtype=object)
    v = c.dot(N)                          # d * vector
    Gm = gram(N, d)
    bvv = int(c.dot(Gm).dot(c)) // (d * d)
    av = v.dot(AL)
    bvav = int(v.dot(G).dot(av)) // (d * d)
    if bvv % 8 or bvav % 4 or not any(c % 2): continue
    # M_v : x with b(x,v), b(x,al v) even  (on basis rows)
    f1 = [int(N[i].dot(G).dot(v)) // (d * d) % 2 for i in range(24)]
    f2 = [int(N[i].dot(G).dot(av)) // (d * d) % 2 for i in range(24)]
    K, k = flint.nmod_mat([f1, f2], 2).nullspace()
    # integer generators of M_v: lifts of kernel vectors + 2*basis
    rows = [[int(K[i, j]) for i in range(24)] for j in range(k)]
    gensMv = [sum(r[i] * N[i] for i in range(24)) for r in rows] + [2 * N[i] for i in range(24)]
    d2 = 2 * d
    newg = [2 * g for g in gensMv] + [v, av]           # everything over denominator 2d
    Hn = flint.fmpz_mat([[int(x) for x in g] for g in newg]).hnf()
    N2 = np.array([[int(x) for x in r] for r in Hn.tolist() if any(int(x) for x in r)], dtype=object)
    if N2.shape[0] != 24: continue
    # remove common factor 2 if possible
    while all(int(x) % 2 == 0 for x in N2.flatten()) and d2 % 2 == 0:
        N2 = N2 // 2; d2 //= 2
    if not valid(N2, d2): continue
    N2 = reduce_basis(N2)
    s2, m2 = score(N2, d2)
    if s2 <= cur:
        N, d, cur, mn = N2, d2, s2, m2; acc += 1
        if cur < best:
            best = cur
            print(f"step {st}: {cur} short pairs (min {mn}), denominator {d}  [{time.time()-t0:.0f}s]", flush=True)
            json.dump(dict(short_pairs=cur, min=mn, d=d, N=[[int(x) for x in r] for r in N]), open(out, "w"))
            if cur == 0:
                print("FOUND extremal J", flush=True); break
print(f"done: best {best}, accepted moves {acc}")
