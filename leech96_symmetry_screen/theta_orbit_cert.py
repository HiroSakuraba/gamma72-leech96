#!/usr/bin/env python3
"""Exact conditional extremal theta coefficients and prime-order orbit filters.

This uses the standard modular-form descriptions for rank-96 even unimodular
lattices and rank-24 even 7-modular lattices. It does not prove existence.
"""
from math import comb
import json,sys

def mul(a,b,n):
    c=[0]*(n+1)
    for i,x in enumerate(a):
        if x:
            for j in range(n-i+1):c[i+j]+=x*b[j]
    return c

def power(a,k,n):
    v=[1]+[0]*n
    for _ in range(k):v=mul(v,a,n)
    return v

def eta_cusp(factors,n):
    v=[1]+[0]*n
    for scale,m in factors:
        for k in range(1,n//scale+1):
            d=scale*k
            t=[0]*(n+1)
            for j in range(min(m,n//d)+1):t[d*j]=(-1)**j*comb(m,j)
            v=mul(v,t,n)
    return [0]+v[:n]

def extremal(basis,n):
    coeff=[1]
    for k in range(1,5):
        assert basis[k][k]==1
        coeff.append(-sum(coeff[j]*basis[j][k] for j in range(k)))
    out=[sum(coeff[j]*basis[j][q] for j in range(5)) for q in range(n+1)]
    assert out[:5]==[1,0,0,0,0]
    return out,coeff

def rank96(n=8):
    e=[1]+[240*sum(d**3 for d in range(1,k+1) if k%d==0)
           for k in range(1,n+1)]
    delta=eta_cusp([(1,24)],n)
    basis=[mul(power(e,12-3*j,n),power(delta,j,n),n)
           for j in range(5)]
    return extremal(basis,n)

def partner7(n=8):
    theta=[0]*(n+1)
    for x in range(-2*n,2*n+1):
        for y in range(-2*n,2*n+1):
            q=x*x+x*y+2*y*y
            if q<=n:theta[q]+=1
    delta=eta_cusp([(1,3),(7,3)],n)
    basis=[mul(power(theta,12-3*j,n),power(delta,j,n),n)
           for j in range(5)]
    return extremal(basis,n)

def main():
    a,ca=rank96()
    b,cb=partner7()
    assert a[5]==565866362880 and a[5]%97==12
    assert b[5]==84672 and b[5]%5==2 and b[5]%13==3
    record={
        "status":"conditional theta series; exact orbit divisibility tests",
        "rank96":{"minimum_if_extremal":10,"norm10_count":a[5],
                  "mod97":a[5]%97,"result":"No automorphism of order 97"},
        "J24_7modular":{"minimum_if_extremal":10,"norm10_count":b[5],
                      "mod5":b[5]%5,"mod13":b[5]%13,
                      "result":"No fixed-point-free order-5 or order-13 automorphism; cyclotomic rank-one Q(zeta_35) and Q(zeta_70) excluded"},
        "basis_coefficients":{"rank96":ca,"J24":cb},
        "coefficients_q0_q8":{"rank96":a,"J24":b},
        "argument":"Any fixed-point-free isometry of prime order p acts on the nonzero norm-10 vectors in orbits of length p.",
        "references":["Quebbemann, Modular lattices in Euclidean spaces, JNT 54 (1995), DOI 10.1006/jnth.1995.1111",
                      "Quebbemann, Atkin-Lehner eigenforms and strongly modular lattices, Enseignement Math. 43 (1997), 55-65"]}
    if len(sys.argv)>1:
        with open(sys.argv[1],"w") as f:json.dump(record,f,indent=2);f.write("\n")
    print(json.dumps(record,indent=2))

if __name__=="__main__":main()
