"""Would-be theta series of an extremal 7-modular lattice of dimension 24 (min 10).
Quebbemann: theta series of even 7-modular lattices of dim 2k lie in C[theta_A, Delta_7] with
theta_A = theta of the binary form [[2,1],[1,4]] (weight 1) and Delta_7 = (eta(z) eta(7z))^3 (weight 3).
Weight 12: f = sum_{i=0}^{4} c_i theta_A^(12-3i) Delta_7^i ; impose f = 1 + 0 q + 0 q^2 + 0 q^3 + 0 q^4 + O(q^5)."""
from fractions import Fraction as Fr
N = 30
def mul(a, b):
    c = [Fr(0)] * N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i): c[i + j] += x * b[j]
    return c
def power(a, k):
    r = [Fr(1)] + [Fr(0)] * (N - 1)
    for _ in range(k): r = mul(r, a)
    return r
th = [Fr(0)] * N
for x in range(-20, 21):
    for y in range(-20, 21):
        n = x * x + x * y + 2 * y * y
        if n < N: th[n] += 1
d = [Fr(0)] * N; d[1] = Fr(1)
for n in range(1, N):
    for _ in range(3):
        d = [d[i] - (d[i - n] if i >= n else 0) for i in range(N)]
    if 7 * n < N:
        for _ in range(3):
            d = [d[i] - (d[i - 7 * n] if i >= 7 * n else 0) for i in range(N)]
basis = [mul(power(th, 12 - 3 * i), power(d, i)) for i in range(5)]
# solve for c: coefficients of q^0..q^4 = (1,0,0,0,0)
import sympy
M = sympy.Matrix([[basis[i][n] for i in range(5)] for n in range(5)])
c = M.solve(sympy.Matrix([1, 0, 0, 0, 0]))
f = [sum(Fr(int(sympy.fraction(c[i])[0]), int(sympy.fraction(c[i])[1])) * basis[i][n] for i in range(5)) for n in range(N)]
print("extremal 7-modular dim-24 theta (q^n counts norm 2n):")
print([str(x) for x in f[:12]])
print("negative:", [n for n in range(N) if f[n] < 0][:5], " non-integral:", [n for n in range(N) if f[n].denominator != 1][:5])
