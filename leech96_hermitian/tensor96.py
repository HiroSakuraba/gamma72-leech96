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

D = sys.argv[1]
leech = json.load(open(f"{D}/leech23.json"))
partners = json.load(open(sys.argv[2]))
which_P = [int(x) for x in sys.argv[3].split(",")]
out_path = sys.argv[4]
max_block = int(sys.argv[5]) if len(sys.argv) > 5 else 30

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

def kill_test(Yamb, p, max_block):
    scale = 8 * p
    A = IntegerMatrix.from_matrix([[int(v) for v in r] for r in Yamb])
    LLL.reduction(A)
    def shortest(A):
        best = None
        for i in range(A.nrows):
            v = [A[i, j] for j in range(A.ncols)]
            n2 = sum(x * x for x in v)
            if best is None or n2 < best[0]:
                best = (n2, v)
        return best
    hist = []
    b = shortest(A)
    hist.append(("LLL", b[0] / scale)); print("    LLL shortest", b[0] / scale, flush=True)
    for bs in [10, 20, 30, 40, 50, 60]:
        if bs > max_block:
            break
        t = time.time()
        par = BKZ.Param(block_size=bs, max_loops=8,
                        flags=BKZ.AUTO_ABORT | BKZ.MAX_LOOPS)
        BKZReduction(A)(par)
        b = shortest(A)
        hist.append((f"BKZ-{bs}", b[0] / scale, round(time.time() - t, 1)))
        print("   ", hist[-1], flush=True)
        if b[0] / scale < 10:
            break
    return b, hist

results = []
for pi in which_P:
    P = partners[pi]
    for si, st in enumerate(leech):
        if "sign" in P and st["sign"] != P["sign"]:
            continue
        t0 = time.time()
        Y, Gm, det, even, Yamb, p = build(P, st)
        print(f"P#{pi} p={P['p']} (trace-min {round(P['trace_min'])}, C1={P['C1']}) x Leech-structure sign {st['sign']}: "
              f"det={det} even={even}", flush=True)
        assert det == 1 and even
        b, hist = kill_test(Yamb, p, max_block)
        found = b[0] / (8 * p)
        status = "KILLED (norm < 10 found)" if found < 10 else "survives this test"
        print(f"  shortest found norm {found}  -> {status}  ({time.time()-t0:.0f}s)", flush=True)
        results.append(dict(P_index=pi, p=P["p"], C1=P["C1"], leech_sign=st["sign"], shortest_norm_found=found,
                            witness_ambient=[int(x) for x in b[1]], history=hist))
        json.dump(results, open(out_path, "w"))
