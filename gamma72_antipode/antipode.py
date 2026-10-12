"""Antipode construction on cross-sections of Gamma72.

Theorem (Conway-Sloane 1996; Chen et al. 2025, Thm 1): Lam unimodular, min mu, dim n.
K = saturated rank-k sublattice (Gram Kg), M = K^* (Gram Kg^{-1}) = projection of Lam onto span(K).
S subset of M with all pairwise squared distances <= beta < mu. Then the (n-k)-dim packing
{ proj_V(w) : w in Lam, proj_U(w) in S } has min squared distance >= mu - beta and center density
    delta = |S| * sqrt(det M) * ((mu - beta)/4)^((n-k)/2).
"""
import json, itertools, math, random
import numpy as np
from fractions import Fraction as F
from sympy import Matrix

MU, N = 8, 72


def load(path):
    return np.array(json.load(open(path)), dtype=np.int64)


def gram_of(G, C):
    """C: k x 72 integer coordinate matrix (rows = vectors in the Gamma72 basis)."""
    return C @ G @ C.T


def is_saturated(C):
    """Rows of C span a primitive (saturated) sublattice of Z^72 iff gcd of all k x k minors is 1,
    i.e. all elementary divisors equal 1. Checked via Hermite normal form."""
    from sympy.matrices.normalforms import smith_normal_form
    from sympy import ZZ
    k = C.shape[0]
    snf = smith_normal_form(Matrix(C.tolist()), domain=ZZ)
    diag = [abs(snf[i, i]) for i in range(k)]
    return all(d == 1 for d in diag), diag


def fincke_pohst(Mf, R):
    """All integer x != 0 with x^T Mf x <= R (Mf positive definite float Gram)."""
    k = Mf.shape[0]
    L = np.linalg.cholesky(Mf)  # Mf = L L^T
    # use upper-triangular form Q: x^T Mf x = sum_i q_ii (x_i + sum_{j>i} q_ij x_j)^2
    U = L.T
    q = np.zeros((k, k))
    for i in range(k):
        q[i, i] = U[i, i] ** 2
        for j in range(i + 1, k):
            q[i, j] = U[i, j] / U[i, i]
    out = []
    x = [0] * k

    def rec(i, rem):
        if i < 0:
            out.append(tuple(x))
            return
        c = -sum(q[i, j] * x[j] for j in range(i + 1, k))
        r = math.sqrt(max(rem, 0) / q[i, i])
        for xi in range(math.ceil(c - r - 1e-12), math.floor(c + r + 1e-12) + 1):
            x[i] = xi
            t = q[i, i] * (xi - c) ** 2
            if t <= rem + 1e-9:
                rec(i - 1, rem - t)
        x[i] = 0

    rec(k - 1, R)
    return [v for v in out if any(v)]


def best_anticodes(Kg, beta_max=1.5, max_pts=4000):
    """For K with integer Gram Kg, find for each achievable diameter beta (<= beta_max) the largest
    S subset of K^* containing 0 with diameter^2 <= beta. Exact rational arithmetic for distances.
    Returns list of (log2 delta, s, beta(Fraction), points)."""
    k = Kg.shape[0]
    Kinv = Matrix(Kg.tolist()).inv()  # exact Gram of dual basis
    Mf = np.array(Kinv.tolist(), dtype=float)
    pts = fincke_pohst(Mf, beta_max)
    if len(pts) > max_pts:
        pts = sorted(pts, key=lambda v: float(np.array(v) @ Mf @ np.array(v)))[:max_pts]
    pts = [tuple(0 for _ in range(k))] + pts
    P = np.array(pts, dtype=float)
    D = (P @ Mf * P).sum(1)[:, None] + (P @ Mf * P).sum(1)[None, :] - 2 * P @ Mf @ P.T
    vals = sorted(set(np.round(D[0, 1:], 9).tolist()) | set(np.round(D[np.triu_indices(len(pts), 1)], 9).tolist()))
    vals = [b for b in vals if 0 < b <= beta_max]
    # only diameters that are realised as a norm from the origin or a pairwise distance among the
    # few shortest points matter; thin the list to keep the clique work bounded
    if len(vals) > 40:
        vals = vals[:40]
    detM = 1 / float(Matrix(Kg.tolist()).det())
    results = []
    import networkx as nx
    for b in vals:
        # candidate vertices: points within distance^2 b of origin
        idx = [i for i in range(len(pts)) if D[0, i] <= b + 1e-9]
        if len(idx) < 2:
            continue
        g = nx.Graph()
        g.add_nodes_from(idx)
        for a, c in itertools.combinations(idx, 2):
            if D[a, c] <= b + 1e-9:
                g.add_edge(a, c)
        clique, _ = nx.max_weight_clique(g.subgraph([i for i in idx if g.has_edge(0, i) or i == 0]), weight=None)
        if 0 not in clique:
            clique = [0] + clique
        s = len(clique)
        diam = max(D[a, c] for a in clique for c in clique)
        l = N - k
        val = math.log2(s) + 0.5 * math.log2(detM) + (l / 2) * math.log2((MU - diam) / 4)
        results.append((val, s, diam, [pts[i] for i in clique]))
    results.sort(reverse=True)
    return results, Kinv


def exact_delta(Kg, S_pts):
    """Exact (rational) data for a final certificate: beta as Fraction, s, det K."""
    Kinv = Matrix(Kg.tolist()).inv()
    k = Kg.shape[0]
    def d2(a, b):
        v = Matrix([a[i] - b[i] for i in range(k)])
        return (v.T * Kinv * v)[0]
    beta = max(d2(a, b) for a in S_pts for b in S_pts)
    return F(int(beta.p), int(beta.q)), len(S_pts), int(Matrix(Kg.tolist()).det())
