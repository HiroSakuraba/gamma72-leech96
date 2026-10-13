#!/usr/bin/env python3
"""Independent pure-Python exact reference for five nontrivial support-six cases.

These are *spot-checks* of the AVX2 algorithm, not a second exhaustive
64-mask enumeration. Every selected support contains at least one residue
pattern passing the modulo-8 bound. The expected (survivors,lifts,short)
counts were computed with direct scalar modular arithmetic.
"""
import itertools
import json
from pathlib import Path
import sys

SQ8=(0,1,4,9,16,9,4,1)
SQ16=(0,1,4,9,16,25,36,49,64,49,36,25,16,9,4,1)
CASES=[
 ((0,2,4,9,12,38),1,(1,6,0)),
 ((0,1,2,3,36,42),47,(1,6,0)),
 ((0,1,6,20,22,34),3,(2,40,0)),
 ((0,1,3,5,21,29),7,(1,9,0)),
 ((0,1,3,12,26,30),15,(1,19,0)),
]

def scalar(A,support,mask):
    k=mask.bit_count()
    minimum_odd_v=48-k if k%2 else k
    opts=[(1,3,5,7) if (mask>>i)&1 else (0,2,4,6) for i in range(6)]
    survivors=lifts=shorts=0
    for r in itertools.product(*opts):
        lower=sum(64 if x==0 else (x*x if i==0 else SQ8[x])
                  for i,x in enumerate(r))
        if lower+minimum_odd_v>128:continue
        cost8=sum(SQ8[sum(r[i]*A[support[i]][j] for i in range(6))%8]
                  for j in range(48))
        if lower+cost8>128:continue
        survivors+=1
        for b in itertools.product((0,1),repeat=5):
            u=[r[0] if r[0] else 8]+[
                (x if x else 8)-(8 if x else 16)*bit
                for x,bit in zip(r[1:],b)]
            u2=sum(x*x for x in u)
            if u2+minimum_odd_v>128:continue
            lifts+=1
            cost16=sum(SQ16[sum(u[i]*A[support[i]][j] for i in range(6))%16]
                       for j in range(48))
            if u2+cost16<=128:shorts+=1
    return survivors,lifts,shorts

if __name__=='__main__':
    if len(sys.argv)!=2:sys.exit(__doc__)
    A=json.loads(Path(sys.argv[1]).read_text())['A']
    for support,mask,expected in CASES:
        actual=scalar(A,support,mask)
        print('mask',mask,'support',support,'survivors/lifts/short',actual,flush=True)
        assert actual==expected,(actual,expected)
    print('Independent reference cases passed')
