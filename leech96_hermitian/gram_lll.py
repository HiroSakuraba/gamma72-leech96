"""Exact LLL on a Gram matrix (small dimension). Returns transform U (rows) with U G U^T reduced."""
from fractions import Fraction as Fr

def gram_lll(G, delta=Fr(99, 100)):
    n = len(G)
    G = [[Fr(x) for x in r] for r in G]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    def gso():
        mu = [[Fr(0)] * n for _ in range(n)]
        Bn = [Fr(0)] * n
        for i in range(n):
            for j in range(i):
                s = G[i][j] - sum(mu[j][k] * mu[i][k] * Bn[k] for k in range(j))
                mu[i][j] = s / Bn[j]
            Bn[i] = G[i][i] - sum(mu[i][k] ** 2 * Bn[k] for k in range(i))
        return mu, Bn
    def addrow(i, j, c):  # b_i += c b_j
        U[i] = [a + c * b for a, b in zip(U[i], U[j])]
        for k in range(n):
            G[i][k] += c * G[j][k]
        for k in range(n):
            G[k][i] = G[i][k] if k != i else G[k][i]
        G[i][i] = G[i][i]  # recompute diag properly below
    def recompute_row(i, old_row_i, j, c):
        pass
    k = 1
    mu, Bn = gso()
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                # b_k -= q b_j ; update Gram exactly
                Gk = G[k][:]
                newrow = [G[k][t] - q * G[j][t] for t in range(n)]
                newdiag = G[k][k] - 2 * q * G[k][j] + q * q * G[j][j]
                newrow[k] = newdiag
                for t in range(n):
                    G[k][t] = newrow[t]; G[t][k] = newrow[t]
                U[k] = [a - q * b for a, b in zip(U[k], U[j])]
                mu, Bn = gso()
        if Bn[k] >= (delta - mu[k][k - 1] ** 2) * Bn[k - 1]:
            k += 1
        else:
            G[k], G[k - 1] = G[k - 1], G[k]
            for r in G:
                r[k], r[k - 1] = r[k - 1], r[k]
            U[k], U[k - 1] = U[k - 1], U[k]
            mu, Bn = gso()
            k = max(k - 1, 1)
    return U, [[int(x) for x in r] for r in G]
