"""Unique Fricke-invariant weight-4 form for Gamma0(11) with q-expansion 1 + 0q + 0q^2 + O(q^3):
would-be theta series of an extremal 11-modular lattice of dimension 8 (min 6).
Basis of M4^+(Gamma0(11)):  E = (E4(z) + 121 E4(11z))/122,  f^2,  f*(E2(z) - 11 E2(11z)),
with f = eta(z)^2 eta(11z)^2 (newform 11a, Fricke eigenvalue -1) - all three Fricke-invariant.
Also numerically re-verifies Fricke invariance of the three basis forms at a sample point."""
from fractions import Fraction as Fr
import cmath, math

N = 40
def sigma(n, k):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)

def E4(scale=1):
    c = [Fr(0)] * N
    c[0] = Fr(1)
    for n in range(1, N):
        if n % scale == 0:
            c[n] += 240 * sigma(n // scale, 3)
    return c

def E2(scale=1):
    c = [Fr(0)] * N
    c[0] = Fr(1)
    for n in range(1, N):
        if n % scale == 0:
            c[n] += -24 * sigma(n // scale, 1)
    return c

def mul(a, b):
    c = [Fr(0)] * N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                c[i + j] += x * b[j]
    return c

def eta_prod():
    # q * prod (1-q^n)^2 (1-q^{11n})^2
    c = [Fr(0)] * N
    c[1] = Fr(1)
    for n in range(1, N):
        for _ in range(2):
            new = c[:]
            for i in range(N - n):
                new[i + n] -= c[i]
            c = new
        if 11 * n < N:
            for _ in range(2):
                new = c[:]
                for i in range(N - 11 * n):
                    new[i + 11 * n] -= c[i]
                c = new
    return c

f = eta_prod()
E = [(a + 121 * b) / 122 for a, b in zip(E4(1), E4(11))]
E2p = [a - 11 * b for a, b in zip(E2(1), E2(11))]
B1, B2 = mul(f, f), mul(f, E2p)
print("f   =", [int(x) for x in f[:12]])
# solve: E + a*B1 + b*B2 has zero q^1, q^2 coefficients
# q^1: E[1] + a*B1[1] + b*B2[1] = 0 ; q^2: E[2] + a*B1[2] + b*B2[2] = 0
a11, a12, r1 = B1[1], B2[1], -E[1]
a21, a22, r2 = B1[2], B2[2], -E[2]
det = a11 * a22 - a12 * a21
a = (r1 * a22 - a12 * r2) / det
b = (a11 * r2 - r1 * a21) / det
theta = [e + a * x + b * y for e, x, y in zip(E, B1, B2)]
print("extremal form coefficients (q^0..q^15):")
for n in range(16):
    print(f"  q^{n}: {theta[n]}  ({float(theta[n]):.4f})")
neg = [n for n in range(N) if theta[n] < 0]
nonint = [n for n in range(N) if theta[n].denominator != 1]
print("negative coefficients at n =", neg[:10])
print("non-integral coefficients at n =", nonint[:10])

# numeric Fricke check at z = 0.1 + 1.2i / sqrt(11)-ish point
def ev(c, z):
    q = cmath.exp(2j * math.pi * z)
    return sum(complex(x) * q ** n for n, x in enumerate(c))
z = complex(0.05, 0.33)
w = -1 / (11 * z)
for name, c in [("E+", E), ("f^2", B1), ("f*E2'", B2)]:
    lhs = ev(c, w) / (11 ** 2) / (z ** 4) * (11 ** 4) / (11 ** 2)  # (sqrt11 z)^-4 f(-1/(11z)) = 11^-2 z^-4 f(w)
    lhs = ev(c, w) * (11 ** -2) * (z ** -4)
    print(f"Fricke check {name}: f|W / f = {lhs / ev(c, z):.6f}")
