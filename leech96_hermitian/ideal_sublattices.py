"""For split primes p of O = Z[alpha] (alpha^2 - alpha + 6 = 0) and each prime ideal P | p,
compute the sublattice P*Leech = (alpha - r) Leech + p Leech, rescale by 1/p (even unimodular, 24-dim)
and find its minimum (Leech type: 4, Niemeier with roots: 2) via SVP."""
import sys, json
import numpy as np
import flint
from fpylll import IntegerMatrix, LLL
sys.path.insert(0, sys.argv[1]); from svp_util import shortest_vector

D = sys.argv[1]
leech = json.load(open(f"{D}/leech23.json"))
primes = [int(x) for x in sys.argv[2].split(",")]
out = {}
for st in leech:
    F = np.array(st["F"], dtype=np.int64); A = np.array(st["A"], dtype=np.int64); H = np.array(st["H"], dtype=np.int64)
    I = np.eye(24, dtype=np.int64)
    for p in primes:
        roots = [r for r in range(p) if (r * r - r + 6) % p == 0]
        for r in roots:
            gens = np.vstack([(A - r * I), p * I])
            Hn = flint.fmpz_mat([[int(v) for v in row] for row in gens]).hnf()
            Y = np.array([[int(v) for v in row] for row in Hn.tolist()], dtype=object)
            Y = Y[[any(v != 0 for v in row) for row in Y]]
            assert Y.shape[0] == 24
            Gm = Y.dot(F.astype(object)).dot(Y.T)
            det = int(flint.fmpz_mat([[int(v) for v in row] for row in Gm]).det())
            amb = Y.dot(H.astype(object))
            M = IntegerMatrix.from_matrix([[int(v) for v in row] for row in amb])
            LLL.reduction(M)
            n2amb, sv = shortest_vector(M)
            n2 = sum(int(x) ** 2 for x in sv) // 8
            key = f"sign{st['sign']}_p{p}_r{r}"
            out[key] = dict(min_over_p=n2 / p, det=det)
            print(f"Leech structure {st['sign']:+d}, p={p}, prime (alpha-{r}, {p}): det={det} (= p^24? {det == p**24}), "
                  f"min = {n2} -> rescaled min {n2/p}", flush=True)
json.dump(out, open(f"{D}/ideal_minima.json", "w"), indent=1)
