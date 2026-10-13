#!/usr/bin/env python3
"""Independent scalar reference checks of three support-five subproblems.

Pure Python; checks exact modulo-8 necessary bounds and the modulo-16
residue norm on every surviving lift. Not a second full 56-billion-class scan.
Usage: python3 support_mod8_reference.py candidate_3cycle.json
"""
import itertools
import json
from pathlib import Path
import sys

SAMPLES={
    (0,1,2,3,4):(0,0,0),
    (0,1,2,12,17):(2,26,0),
    (43,44,45,46,47):(0,0,0),
}

def reference(A,support):
    sq8=(0,1,4,9,16,9,4,1)
    passed=lifts=short=0
    for r in itertools.product(range(8),repeat=5):
        lower=(64 if r[0]==0 else r[0]*r[0])+sum(
            64 if t==0 else sq8[t] for t in r[1:])
        if lower>128:continue
        for j in range(48):
            v=sum(r[i]*A[support[i]][j] for i in range(5))%8
            lower+=sq8[v]
            if lower>128:break
        if lower>128:continue
        passed+=1
        for bits in itertools.product((0,1),repeat=4):
            u=[r[0] if r[0] else 8]
            for i in range(4):
                t=r[i+1]
                u.append((t if t else 8)-(8 if t else 16)*bits[i])
            u2=sum(x*x for x in u)
            odd=sum(x%2 for x in u)
            if u2+(48-odd if odd%2 else odd)>128:continue
            lifts+=1
            for j in range(48):
                v=sum(u[i]*A[support[i]][j] for i in range(5))%16
                u2+=min(v,16-v)**2
                if u2>128:break
            if u2<=128:short+=1
    return passed,lifts,short

def main():
    if len(sys.argv)!=2:sys.exit(__doc__)
    A=json.loads(Path(sys.argv[1]).read_text())["A"]
    for support,expected in SAMPLES.items():
        got=reference(A,support)
        print("support:",support,"passes/lifts/short:",got,flush=True)
        assert got==expected,(got,expected)
    print("Independent scalar reference checks passed")

if __name__=="__main__":
    main()
