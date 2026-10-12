"""Maximally symmetric Paley Construction-A lattices in dimension 96.
C = order-48 skew conference matrix (Paley, q = 47). A = m C + c I (mod N), centred.
L = (1/sqrt N) {(u, v) : v = u A mod N},  even unimodular iff  A A^T = -I (mod N) and every row of [I|A]
has squared length = 0 (mod 2N).  Any automorphism of C (contains PSL(2,47) acting on both halves) is an
automorphism of L.  min(L) <= N because of (0, N e_i).  Kill test by randomized BKZ."""
import sys, json, time, random
import numpy as np
import flint
from fpylll import IntegerMatrix, LLL, BKZ
from fpylll.algorithms.bkz2 import BKZReduction

q = 47
QR = {(x * x) % q for x in range(1, q)}
chi = lambda a: 0 if a % q == 0 else (1 if a % q in QR else -1)
Q = np.array([[chi(j - i) for j in range(q)] for i in range(q)], dtype=np.int64)
C = np.zeros((48, 48), dtype=np.int64)
C[0, 1:] = 1; C[1:, 0] = -1; C[1:, 1:] = Q
I = np.eye(48, dtype=np.int64)
assert (C.T == -C).all() and (C @ C.T == 47 * I).all()

def centred(x, N):
    r = x % N
    return r - N if r > N // 2 else r

Ns = [int(x) for x in sys.argv[1].split(",")]
top = int(sys.argv[2]); runs = int(sys.argv[3]); out = sys.argv[4]
cands = []
for N in Ns:
    seen = set()
    for m in range(N):
        for c in range(N):
            if (47 * m * m + c * c + 1) % N:
                continue
            mt, ct = centred(m, N), centred(c, N)
            key = (abs(mt), abs(ct))
            if key in seen:
                continue
            A = mt * C + ct * I
            if ((1 + (A * A).sum(1)) % (2 * N) != 0).any():
                continue
            seen.add(key)
            cands.append((N, mt, ct, A))
print("candidates (N, m, c):", [(N, m, c) for N, m, c, _ in cands], flush=True)

res = []
for N, m, c, A in cands:
    G = np.block([[(I + A @ A.T) // N, A], [A.T, N * I]])
    assert ((I + A @ A.T) % N == 0).all()
    det = int(flint.fmpz_mat(G.tolist()).det()); even = (np.diag(G) % 2 == 0).all()
    assert det == 1 and even
    basis = np.vstack([np.hstack([I, A]), np.hstack([0 * I, N * I])])
    best = None
    t0 = time.time()
    for r in range(runs):
        rng = random.Random(31 * N + 7 * r + m)
        M = IntegerMatrix.from_matrix(basis.tolist())
        for _ in range(400):
            i, j = rng.sample(range(96), 2); M[i].addmul(M[j], rng.choice([-1, 1]))
        LLL.reduction(M)
        for bs in [20, 30, 40]:
            if bs > top: break
            BKZReduction(M)(BKZ.Param(block_size=bs, max_loops=4, flags=BKZ.AUTO_ABORT | BKZ.MAX_LOOPS))
            k = min(range(96), key=lambda i: sum(M[i, j] ** 2 for j in range(96)))
            vec = [M[k, j] for j in range(96)]
            nm = sum(x * x for x in vec) / N
            if best is None or nm < best[0]:
                best = (nm, vec)
            if nm < 10: break
        if best[0] < 10: break
    u, v = np.array(best[1][:48]), np.array(best[1][48:])
    member = ((v - u @ A) % N == 0).all()
    print(f"N={N} m={m} c={c}: even unimodular; shortest found norm {best[0]} (member {member}) "
          f"{'KILLED' if best[0] < 10 else 'survives'}  [{time.time()-t0:.0f}s]", flush=True)
    res.append(dict(N=N, m=m, c=c, shortest=best[0], vector=[int(x) for x in best[1]], member=bool(member)))
    json.dump(res, open(out, "w"))
