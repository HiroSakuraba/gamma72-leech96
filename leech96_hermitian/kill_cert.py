"""Explicit kill certificate via the ideal-class obstruction.
For candidate L = P (x)_O Leech (O = Z[alpha], alpha^2-alpha+6=0, class number 3):
  q = prime ideal (alpha - r, 2) whose twist qLeech has rescaled minimum 2 (a 'bad' class),
  x = shortest vector of qLeech (Leech norm 4),
  v = shortest vector of conj(q) P in the form h(v,v)/4  (so v/2 runs over q^{-1}P),
  z = (v/2) (x) x  = blocks ( (v_i * x) / 2 )_i  in Leech^4.
Check z in L exactly and report its norm  (= h_P(v/2) * |x|^2 <= 8)."""
import sys, json
import numpy as np
import flint
sys.path.insert(0, sys.argv[1])
from svp_util import shortest_vector
from fpylll import IntegerMatrix
import importlib.util
spec = importlib.util.spec_from_file_location("t96", sys.argv[1] + "/tensor96_lib.py")
t96 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t96)
from partners import fincke_pohst_min, T2

D = sys.argv[1]
cands = json.load(open(sys.argv[2]))
leech = {s["sign"]: s for s in json.load(open(f"{D}/leech23.json"))}
imin = json.load(open(f"{D}/ideal_minima.json"))

def omul(s, t):  # (a+b al)(c+d al), al^2 = al - 6
    a, b = s; c, d = t
    return (a * c - 6 * b * d, a * d + b * c + b * d)

def conj(s):
    a, b = s
    return (a + b, -b)

out = []
for pi, P in enumerate(cands):
    st = leech[P["sign"]]
    F = np.array(st["F"], dtype=np.int64); A = np.array(st["A"], dtype=np.int64); H = np.array(st["H"], dtype=np.int64)
    I = np.eye(24, dtype=np.int64)
    r = [rr for rr in (0, 1) if imin[f"sign{P['sign']}_p2_r{rr}"]["min_over_p"] == 2.0][0]
    # q = (alpha - r, 2);  qLeech = (A - rI) Leech + 2 Leech
    gens = np.vstack([A - r * I, 2 * I])
    Hn = flint.fmpz_mat(gens.tolist()).hnf()
    Yq = np.array([[int(v) for v in row] for row in Hn.tolist()], dtype=np.int64)
    Yq = Yq[np.abs(Yq).sum(1) > 0]
    M = IntegerMatrix.from_matrix((Yq @ H).tolist())
    n2, xamb = shortest_vector(M)
    x_leech = np.array(flint.fmpz_mat([list(map(int, xamb))]).tolist()[0])  # ambient
    # Leech coordinates of x: solve x = c H
    c = flint.fmpq_mat([list(map(int, xamb))]) * flint.fmpq_mat(H.tolist()).inv()
    xc = np.array([int(str(v)) for v in c.entries()], dtype=np.int64)
    assert xc @ F @ xc == n2 // 8 == 4
    # conj(q) = (conj(alpha - r), 2) = (1 - alpha - r, 2)
    qbar_gens = [(1 - r, -1), (2, 0)]
    Bp = np.array(P["B"], dtype=np.int64)   # rows (a1,b1,...,a4,b4)
    rows = []
    for row in Bp:
        for g in qbar_gens:
            rows.append([x for i in range(4) for x in omul(g, (row[2 * i], row[2 * i + 1]))])
    Hn2 = flint.fmpz_mat([[int(t) for t in rr] for rr in rows]).hnf()
    Bq = np.array([[int(v) for v in rr] for rr in Hn2.tolist()], dtype=np.int64)
    Bq = Bq[np.abs(Bq).sum(1) > 0]
    T8 = np.kron(np.eye(4, dtype=np.int64), T2)
    Gq = Bq @ T8 @ Bq.T          # trace form of sum v_i conj(w_i)  (before /p and /4)
    # Hermitian norm of v/2 in P's form: h = (1/p) sum |v_i/2|^2 = trace/(2*4*p)
    # shortest v by enumeration on Gq
    pass
    # simple Fincke-Pohst in 8 dims with the vector returned
    import math
    best = None
    Gf = Gq.astype(float)
    L = np.linalg.cholesky(Gf); U = L.T
    k = 8; q = np.zeros((k, k))
    for i in range(k):
        q[i, i] = U[i, i] ** 2
        for j in range(i + 1, k):
            q[i, j] = U[i, j] / U[i, i]
    bound = 2 * 4 * P["p"] * 4.0 + 1e-6   # look for h(v/2) <= 4
    sol = []
    xx = [0] * k
    def rec(i, rem):
        if i < 0:
            if any(xx): sol.append(list(xx))
            return
        cc = -sum(q[i, j] * xx[j] for j in range(i + 1, k))
        rr = math.sqrt(max(rem, 0) / q[i, i])
        for xi in range(math.ceil(cc - rr - 1e-9), math.floor(cc + rr + 1e-9) + 1):
            xx[i] = xi
            t = q[i, i] * (xi - cc) ** 2
            if t <= rem + 1e-9: rec(i - 1, rem - t)
        xx[i] = 0
    rec(k - 1, bound)
    coeffs = min(sol, key=lambda s: np.array(s) @ Gq @ np.array(s))
    v = np.array(coeffs) @ Bq     # (a1,b1,...,a4,b4) in conj(q)P
    hv2 = (v @ T8 @ v) / (2 * 4 * P["p"])
    # z blocks: (v_i * x)/2 with s*x = x M_s, M_s = a I + b A  (Leech coords)
    blocks = []
    for i in range(4):
        a, b = int(v[2 * i]), int(v[2 * i + 1])
        w = xc @ (a * I + b * A)
        assert (w % 2 == 0).all(), "v_i x not in 2 Leech"
        blocks.append(w // 2)
    z = np.concatenate(blocks)
    norm = sum(int(blocks[i] @ F @ blocks[i]) for i in range(4)) / P["p"]
    # membership in L
    Y, Gm, det, even, Yamb, p = t96.build(P, st)
    sol_c = flint.fmpq_mat([list(map(int, z))]) * flint.fmpq_mat([[int(t) for t in row] for row in Y]).inv()
    member = all(int(e.q) == 1 for e in sol_c.entries())
    print(f"P#{pi} p={P['p']} sign {P['sign']:+d}: bad ideal (alpha-{r},2); |x|^2={n2//8}; "
          f"h(v/2)={hv2}; z in L: {member}; norm(z) = {norm}", flush=True)
    out.append(dict(P_index=pi, p=P["p"], sign=P["sign"], ideal=f"(alpha-{r},2)", h_v=hv2, x_norm=4,
                    in_L=member, norm=norm, z_leech_coords=[int(t) for t in z]))
json.dump(out, open(sys.argv[3], "w"))
print("all killed:", all(o["in_L"] and o["norm"] < 10 for o in out))
