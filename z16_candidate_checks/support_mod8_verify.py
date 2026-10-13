#!/usr/bin/env python3
"""Exact finite radius-8 scans for the published 96-dimensional Z/16 candidate.

Usage: python3 support_mod8_verify.py candidate_3cycle.json five|six result.json
five: EVERY vector with exactly five nonzero centered first-half u coordinates.
six: every vector with exactly six nonzero centered u coordinates, all odd.
Requires a C++17 compiler with OpenMP and an x86-64 AVX2 CPU.
Neither test establishes the full 96-dimensional minimum norm.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import tempfile
from math import comb

EXPECTED = {
    "five": {"support":5,"supports":1712304,"mod8_classes":56108777472,
             "pass_mod8":15278,"lifts":127019,
             "passes_by_odd":[0,5714,635,3895,1498,3536],
             "lifts_by_odd":[0,47344,8845,30077,18427,22326]},
    "six": {"support":6,"supports":12271512,"mod8_classes":50264113152,
            "pass_mod8":14898,"lifts":252059,
            "passes_by_odd":[0,0,0,0,0,0,14898],
            "lifts_by_odd":[0,0,0,0,0,0,252059]},
}

def validate(A, sha):
    assert len(A)==48 and all(len(row)==48 for row in A)
    assert all(type(x)==int and -8<=x<=8 for row in A for x in row)
    assert all(A[i][j]%2==int(i!=j) for i in range(48) for j in range(48))
    for i in range(48):
        assert (1+sum(x*x for x in A[i]))%32==0
        for j in range(i+1):
            assert (sum(A[i][k]*A[j][k] for k in range(48)) + int(i==j))%16==0
    digest=hashlib.sha256(b"".join(int(x).to_bytes(8,"little",signed=True)
                                     for row in A for x in row)).hexdigest()
    assert sha==digest,"Candidate matrix checksum mismatch"
    return digest

def main():
    if len(sys.argv)!=4 or sys.argv[2] not in EXPECTED:sys.exit(__doc__)
    C=json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    A=C["A"]; sha=validate(A,C["A_sha256"])
    scope=sys.argv[2]; ex=EXPECTED[scope]; k=ex["support"]
    if platform.machine() not in ("x86_64","AMD64"):
        raise RuntimeError("This optimized verifier requires an AVX2-capable x86-64 CPU")
    code=Path(__file__).with_name(
        "support_mod8_exact.cpp" if scope=="five" else "support6_allodd.cpp")
    with tempfile.TemporaryDirectory(prefix="leech96_mod8_") as tmp:
        tmp=Path(tmp)
        header='#include <cstdint>\nstatic const int8_t A[48][48]={\n'
        header+=',\n'.join('{'+','.join(map(str,row))+'}' for row in A)
        header+='\n};\n'
        (tmp/"matrix48.h").write_text(header,encoding="utf-8")
        exe=tmp/"scan"
        cmd=[os.environ.get("CXX","g++"),"-O3","-std=c++17","-mavx2",
             "-fopenmp","-I",str(tmp)]
        if scope=="five":cmd.append("-DSUPPORT=5")
        cmd.extend([str(code),"-o",str(exe)])
        subprocess.run(cmd,check=True)
        env={**os.environ,"OMP_NUM_THREADS":os.environ.get("OMP_NUM_THREADS","4")}
        run=subprocess.run([str(exe),"0",str(49-k)],text=True,
                            capture_output=True,env=env)
    print(run.stdout,flush=True)
    pattern=(r"SUMMARY (?:support=(\d+) )?first=\[(\d+),(\d+)\) "
             r"supports=(\d+) mod8_classes=(\d+) pass_mod8=(\d+) "
             r"lifts=(\d+) shorts=(\d+) best_num=(\d+)")
    m=re.search(pattern,run.stdout)
    assert m,(run.returncode,run.stderr,run.stdout)
    mode,start,end,supports,classes,passed,lifts,short,best=m.groups()
    assert mode is None or int(mode)==k
    assert int(start)==0 and int(end)==49-k
    assert int(supports)==comb(48,k)==ex["supports"]
    assert (int(classes),int(passed),int(lifts))==(
        ex["mod8_classes"],ex["pass_mod8"],ex["lifts"])
    assert int(short)==0 and run.returncode==0,"Short vector found: construction eliminated"
    pass_buckets=ex["passes_by_odd"];lift_buckets=ex["lifts_by_odd"]
    if scope=="five":
        mm=re.search(r"MOD8_PASSES_BY_ODD_COUNTS((?: \d+)+)",run.stdout)
        nn=re.search(r"LIFTS_BY_ODD_COUNTS((?: \d+)+)",run.stdout)
        assert mm and nn
        assert list(map(int,mm.group(1).split()))==pass_buckets
        assert list(map(int,nn.group(1).split()))==lift_buckets
    result={
        "candidate_A_sha256":sha,
        "scope":("exactly 5 nonzero centered first-half coordinates, any parity"
                 if scope=="five" else "exactly 6 nonzero centered first-half coordinates, all odd"),
        "norm_squared_threshold":8,
        "exact_support_count":int(supports),
        "mod8_classes_evaluated":int(classes),
        "mod8_survivors":int(passed),
        "admissible_mod16_lifts_exactly_checked":int(lifts),
        "mod8_survivors_by_odd_coordinate_count":pass_buckets,
        "admissible_lifts_by_odd_coordinate_count":lift_buckets,
        "short_vectors_found":int(short),
        "complete_within_stated_scope":True,
        "full_96_dimensional_minimum_certified":False}
    Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
