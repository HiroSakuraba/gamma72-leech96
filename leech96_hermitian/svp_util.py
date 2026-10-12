"""Exact shortest vector via fpylll enumeration (no strategy files needed)."""
from fpylll import IntegerMatrix, LLL, GSO, BKZ
from fpylll.fplll.enumeration import Enumeration, EnumerationError
from fpylll.algorithms.bkz2 import BKZReduction


def reduce(M, block=20):
    LLL.reduction(M)
    if block and M.nrows > block:
        BKZReduction(M)(BKZ.Param(block_size=block, max_loops=4, flags=BKZ.AUTO_ABORT | BKZ.MAX_LOOPS))
    return M


def shortest_vector(M, block=20):
    """Return (squared norm, coefficient-free ambient vector) of a shortest nonzero vector, exact."""
    reduce(M, block)
    G = GSO.Mat(M)
    G.update_gso()
    n = M.nrows
    R = min(sum(M[i, j] ** 2 for j in range(M.ncols)) for i in range(n))
    enum = Enumeration(G)
    sols = enum.enumerate(0, n, R * 1.0001, 0)
    best = None
    for dist, coeffs in sols:
        v = [sum(int(round(coeffs[i])) * M[i, j] for i in range(n)) for j in range(M.ncols)]
        n2 = sum(x * x for x in v)
        if n2 > 0 and (best is None or n2 < best[0]):
            best = (n2, v)
    return best


def count_short(M, bound, block=20, max_sols=100000):
    """Enumerate vectors with squared ambient norm <= bound (up to sign); returns list of vectors."""
    reduce(M, block)
    G = GSO.Mat(M)
    G.update_gso()
    n = M.nrows
    enum = Enumeration(G, nr_solutions=max_sols)
    try:
        sols = enum.enumerate(0, n, bound * 1.0001, 0)
    except EnumerationError:
        return []
    out = []
    for dist, coeffs in sols:
        v = [sum(int(round(coeffs[i])) * M[i, j] for i in range(n)) for j in range(M.ncols)]
        if any(v):
            out.append(v)
    return out
