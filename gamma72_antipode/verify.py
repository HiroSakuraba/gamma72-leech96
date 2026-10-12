"""Independent verification of the antipode certificates (does not import the search code).

For each certificate (C, S):
  1. K = C G C^T matches the recorded Gram; rows of C are a saturated set in Z^72 (Smith form all ones).
  2. For each s in S, solve C G w = s for an INTEGER w (i.e. exhibit w in Gamma72 whose projection
     onto U = span(K) is the dual-lattice point with coordinates s). Existence of all these w is what
     makes the |S| translates actually lie in Gamma72.
  3. beta = max_{s,t in S} (s-t)^T K^{-1} (s-t) computed exactly; delta recomputed exactly.
  4. Spot-check of the distance bound on short vectors: for every vector x of norm <= 12 in a pool of
     Gamma72 vectors whose projection lies in S - S, check |proj_V x|^2 >= 8 - beta.
The packing bound itself is the elementary inequality |pi_V x|^2 = |x|^2 - |pi_U x|^2 >= 8 - beta
for x in Gamma72 \\ 0 with pi_U x in S - S, using min(Gamma72) = 8 (Nebe 2012)."""
import sys, json
from fractions import Fraction as F
import numpy as np
from sympy import Matrix, ZZ, Rational, sqrt, log, N as Nf
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form

G = np.array(json.load(open(sys.argv[1])), dtype=np.int64)
certs = json.load(open(sys.argv[2]))
GM = Matrix(G.tolist())
assert GM.det() == 1

def solve_integer(A, b):
    """Integer solution of A w = b (A: k x n full row rank, integer). Uses Smith form A = U^-1 D V^-1."""
    from sympy.matrices.normalforms import smith_normal_decomp
    D, U, V = smith_normal_decomp(A, domain=ZZ)  # D = U A V
    k = A.rows
    c = U * b
    y = Matrix.zeros(A.cols, 1)
    for i in range(k):
        assert c[i] % D[i, i] == 0
        y[i] = c[i] // D[i, i]
    w = V * y
    assert A * w == b
    return w

# short-vector pool for the spot check
n = 72
E = np.eye(n, dtype=np.int64)
pool = [E[i] for i in range(n)]
for i in range(n):
    for j in range(i + 1, n):
        for s in (1, -1):
            v = E[i] + s * E[j]
            if v @ G @ v <= 12:
                pool.append(v)
pool = np.array(pool)
pool = np.concatenate([pool, -pool])

summary = []
for key in sorted(certs, key=int):
    c = certs[key]
    C = Matrix(c["C"]); k = C.rows; l = n - k
    K = C * GM * C.T
    assert K == Matrix(c["K"]), "Gram mismatch"
    snf = smith_normal_form(C, domain=ZZ)
    assert all(abs(snf[i, i]) == 1 for i in range(k)), "not saturated"
    detK = K.det()
    Kinv = K.inv()
    S = [Matrix(p) for p in c["S"]]
    A = C * GM
    ws = [solve_integer(A, s) for s in S]          # lattice vectors realising each translate
    beta = max(((s - t).T * Kinv * (s - t))[0] for s in S for t in S)
    assert beta < 8
    s_count = len(S)
    # exact center density: s * detK^{-1/2} * ((8-beta)/4)^{l/2}
    delta = s_count * detK ** Rational(-1, 2) * ((8 - beta) / 4) ** Rational(l, 2)
    l2 = float(Nf(log(delta, 2), 20))
    # spot check on short vectors
    Sset = {tuple(int(v) for v in s) for s in S}
    diffs = {tuple(a - b for a, b in zip(s, t)) for s in Sset for t in Sset}
    P = pool @ G @ np.array(c["C"]).T  # dual coordinates <x, c_i>
    Kinv_f = np.array(Kinv.tolist(), dtype=float)
    checked = 0
    worst = 1e9
    for x, px in zip(pool, P):
        if tuple(int(v) for v in px) in diffs:
            nx2 = float(x @ G @ x)
            pu2 = float(px @ Kinv_f @ px)
            worst = min(worst, nx2 - pu2)
            checked += 1
    assert worst >= float(8 - beta) - 1e-9
    hermite = {1: 1, 2: Rational(4, 3), 3: 2, 4: 4, 5: 8, 6: Rational(64, 3), 7: 64}[k]
    bare_bound = 2 ** Rational(l, 2) * (hermite / 8 ** k) ** Rational(1, 2)
    l2_bare = float(Nf(log(bare_bound, 2), 20))
    print(f"dim {l}: K saturated, det K = {detK}, |S| = {s_count}, beta = {beta}, "
          f"all {s_count} translates realised by integer vectors of Gamma72")
    print(f"   delta = {s_count} * {detK}^(-1/2) * ({(8-beta)/4})^({l}/2) = {float(delta):.6g} = 2^{l2:.5f}")
    print(f"   best possible bare cross-section of Gamma72 (Hermite bound): 2^{l2_bare:.5f};  "
          f"ratio = {2**(l2-l2_bare):.4f};  spot-check: {checked} short vectors, min |pi_V x|^2 = {worst:.4f} >= {float(8-beta):.4f}")
    summary.append(dict(dim=l, k=k, detK=int(detK), s=s_count, beta=str(beta), delta=float(delta), log2delta=l2,
                        bare_bound_log2=l2_bare, K=c["K"], C=c["C"], S=c["S"],
                        witnesses=[[int(v) for v in w] for w in ws]))
json.dump(summary, open(sys.argv[3], "w"), indent=1)
print("ALL CHECKS PASSED")
