#!/usr/bin/env python3
"""Exact rank-96 lattice support-six certificate for one parity mask.

Usage: python3 support6_mask_verify.py candidate_3cycle.json MASK output.json
Run MASK=0,...,63 for all combinations of odd/even support coefficients.
Requires Python 3.9+, a C++17/OpenMP compiler, and an x86-64 AVX2 CPU.
Each invocation validates the exact 48x48 candidate matrix fingerprint and
Type-II lattice congruences, exhaustively enumerates the given mask, checks
all centered mod-16 survivors, and outputs a machine-readable certificate.
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

KNOWN = {
    # Previously completed independent five-odd / one-even exact scans.
    62:(10290,79364), 61:(9609,76698), 59:(9043,71622),
    55:(9429,74884), 47:(8898,70919), 31:(9153,72346),
    # Previously completed independent one-odd / five-even exact scans.
    1:(503,3018), 2:(557,3682), 4:(650,4339),
    8:(667,4438), 16:(617,4014), 32:(590,3924),
    # Independently published previous six-odd check.
    63:(14898,252059),
}

def check_matrix(c):
    A=c['A']
    assert len(A)==48 and all(len(row)==48 for row in A)
    assert all(type(x) is int and -8<=x<=8 for row in A for x in row)
    sha=hashlib.sha256(b''.join(int(x).to_bytes(8,'little',signed=True)
                            for row in A for x in row)).hexdigest()
    assert sha==c['A_sha256'],(sha,c['A_sha256'])
    for i in range(48):
        assert (1+sum(x*x for x in A[i]))%32==0
        for j in range(i+1):
            assert (sum(A[i][k]*A[j][k] for k in range(48))+(i==j))%16==0
        for j in range(48):
            assert A[i][j]%2==int(i!=j)
    return A,sha

def main():
    if len(sys.argv)!=4:sys.exit(__doc__)
    candidate=json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))
    A,sha=check_matrix(candidate)
    mask=int(sys.argv[2])
    if not 0<=mask<=63:raise ValueError("MASK must be in 0..63")
    if platform.machine().lower() not in ('x86_64','amd64'):
        raise RuntimeError("Optimized exact verifier requires x86-64 AVX2")
    source=Path(__file__).with_name('support6_mask_exact.cpp')
    with tempfile.TemporaryDirectory(prefix='leech96_support6_') as t:
        tmp=Path(t)
        header=('#include <cstdint>\nstatic const int8_t A[48][48]={\n'
               +',\n'.join('{'+','.join(map(str,row))+'}' for row in A)
               +'\n};\n')
        (tmp/'matrix48.h').write_text(header)
        exe=tmp/'scan'
        subprocess.run([os.environ.get('CXX','g++'),'-std=c++17','-O3',
                        '-mavx2','-fopenmp','-I',str(tmp),str(source),
                        '-o',str(exe)],check=True)
        process=subprocess.run([str(exe),str(mask),'0','43'],capture_output=True,
            text=True,env={**os.environ,'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS','4')})
    print(process.stdout,flush=True)
    match=re.search(r'SUMMARY oddmask=(\d+) first=\[(\d+),(\d+)\) supports=(\d+) classes=(\d+) survivors=(\d+) lifts=(\d+) short=(\d+) best_num=(\d+)',process.stdout)
    assert match,(process.returncode,process.stdout,process.stderr)
    m,start,end,supports,classes,survive,lifts,short,best=map(int,match.groups())
    assert (m,start,end)==(mask,0,43)
    assert supports==comb(48,6)
    assert classes==supports*4096
    if mask in KNOWN:
        assert (survive,lifts)==KNOWN[mask],(mask,(survive,lifts),KNOWN[mask])
    if short:
        assert process.returncode==3
    else:
        assert process.returncode==0,(process.returncode,process.stderr)
    result={
        'candidate_A_sha256':sha,
        'scope':'first-half vector with exactly six nonzero centered coordinates; parity in ascending support order given by mask',
        'oddmask':mask,
        'number_of_odd_support_coordinates':mask.bit_count(),
        'supports':supports,
        'residue_patterns_mod8':classes,
        'mod8_survivors':survive,
        'exact_mod16_lifts_checked':lifts,
        'short_vectors_norm_sq_at_most_8':short,
        'best_witness_numerator':None if best==2147483647 else best,
        'complete_for_mask':True,
        'global_minimum_certified':False
    }
    Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if short:sys.exit(3)

if __name__=='__main__':main()
