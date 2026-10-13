#!/usr/bin/env python3
"""Portable exact scalar references (spots only, not exhaustive 1.2T scan).
Check examples with positive and zero modulo-8 survivors for all-odd support 7.
"""
import itertools,json,sys
from pathlib import Path
S8=(0,1,4,9,16,9,4,1)
S16=(0,1,4,9,16,25,36,49,64,49,36,25,16,9,4,1)
CASES={
 (39,41,43,44,45,46,47):(1,7,0),
 (38,39,40,41,42,43,44):(0,0,0),
}
def check(A,ind):
    p8=lifts=short=0
    for r in itertools.product((1,3,5,7),repeat=7):
        umin=r[0]*r[0]+sum(S8[t] for t in r[1:])
        if umin+41>128:continue
        residues=[sum(r[i]*A[ind[i]][j] for i in range(7))%8 for j in range(48)]
        if umin+sum(S8[t] for t in residues)>128:continue
        p8+=1
        for bits in itertools.product((0,1),repeat=6):
            u=[r[0]]+[r[i]-8*bits[i-1] for i in range(1,7)]
            n=sum(t*t for t in u)
            if n+41>128:continue
            lifts+=1
            n+=sum(S16[sum(u[i]*A[ind[i]][j] for i in range(7))%16]
                   for j in range(48))
            if n<=128:short+=1
    return p8,lifts,short
if __name__=="__main__":
    if len(sys.argv)!=2:sys.exit(__doc__)
    A=json.loads(Path(sys.argv[1]).read_text())["A"]
    for indices,expected in CASES.items():
        got=check(A,indices)
        print(indices,got,flush=True)
        assert got==expected
    print("Scalar reference cases passed")
