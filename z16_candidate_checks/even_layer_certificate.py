#!/usr/bin/env python3
"""Exact, reproducible even-u exclusion for the published Z/16 lattice candidate.

Requires Python 3.9+ and a C++17 compiler. No third-party Python packages.
Does not exclude odd-u vectors. The existing layer4.py handles u=0 mod4.

Usage: python3 even_layer_certificate.py candidate_3cycle.json report.json
"""
import hashlib, json, os, re, subprocess, sys, tempfile
from pathlib import Path
N=48
FULL=(1<<N)-1

def basis(rows):
    piv={}
    for row in rows:
        v=sum(int(b)<<j for j,b in enumerate(row))
        while v:
            p=v.bit_length()-1
            if p in piv:v^=piv[p]
            else:piv[p]=v;break
    return list(piv.values())

def weight_enumerator(rows,offset=0):
    b=basis(rows)
    assert len(b)==24
    counts=[0]*(N+1)
    x=offset
    counts[x.bit_count()]+=1
    for n in range(1,1<<24):
        x^=b[(n&-n).bit_length()-1]
        counts[x.bit_count()]+=1
    assert sum(counts)==1<<24
    return counts

def validate(A):
    assert len(A)==48 and all(len(r)==48 for r in A)
    assert all(type(v) is int and -128<=v<=127 for r in A for v in r)
    assert all(A[i][i]%4==0 for i in range(N))
    assert all((A[i][j]+A[j][i])%4==0 for i in range(N) for j in range(N))
    assert all(A[i][j]%2==int(i!=j) for i in range(N) for j in range(N))
    for i in range(N):
        assert (1+sum(v*v for v in A[i]))%32==0
        for j in range(i+1):
            assert (sum(A[i][k]*A[j][k] for k in range(N))+int(i==j))%16==0

def header(A,M):
    rows=[sum(1<<j for j,v in enumerate(row) if v) for row in M]
    cols=[sum(1<<i for i in range(N) if M[i][j]) for j in range(N)]
    lines=["#include <cstdint>","static const int8_t A[48][48]={"]
    lines += ["{"+",".join(map(str,row))+"}," for row in A]
    lines += ["};",
        "static const uint64_t row_masks[48]={"+",".join(str(v)+"ULL" for v in rows)+"};",
        "static const uint64_t column_masks[48]={"+",".join(str(v)+"ULL" for v in cols)+"};"]
    return "\n".join(lines)+"\n"

def main():
    if len(sys.argv)!=3:sys.exit(__doc__)
    candidate=json.loads(Path(sys.argv[1]).read_text())
    A=candidate["A"]
    validate(A)
    M=[[int(A[i][j]%4==3) for j in range(N)] for i in range(N)]
    Mt=list(map(list,zip(*M)))
    assert all((M[i][j]^M[j][i])==int(i!=j) for i in range(N) for j in range(N))
    assert all(sum(M[i][j] for i in range(N))%2==0 for j in range(N))
    assert all(sum(M[i][j] for j in range(N))%2==1 for i in range(N))
    C0=weight_enumerator(M)
    C1=weight_enumerator(Mt)
    C2=weight_enumerator(M,offset=FULL)
    assert [(c[12],c[14],c[16]) for c in [C0,C1,C2]]==[
        (12972,0,356730),(17296,0,535095),(4324,0,178365)]
    assert C0[11]==4324 and C1[11]==0 and C2[0]==0 and C2[1]==1
    assert all(not c[2*j] for c in (C0,C1,C2) for j in range(1,6))
    assert all(c[12]>0 for c in (C0,C1,C2))
    cpp=Path(__file__).with_name("even_layer_core.cpp")
    with tempfile.TemporaryDirectory(prefix="leech96_even_") as tmp:
        tmp=Path(tmp)
        (tmp/"matrix48.h").write_text(header(A,M))
        exe=tmp/"scan"
        subprocess.run([os.environ.get("CXX","g++"),"-O3","-std=c++17",
            str(cpp),"-I",str(tmp),"-o",str(exe)],check=True)
        result=subprocess.run([str(exe)],text=True,capture_output=True,check=True)
    matches=re.findall(r"family=(\d+) k=(\d+) supports=(\d+) cases=(\d+) modular_pass=(\d+) short=(\d+)",result.stdout)
    assert len(matches)==6,result.stdout
    expected_pass={(0,12):0,(0,16):7,(1,12):0,(1,16):16,(2,12):0,(2,16):3}
    counts=[C0,C1,C2]
    rows=[]
    for fields in matches:
        f,k,s,v,p,w=map(int,fields)
        assert s==counts[f][k]
        assert v==s*(97 if k==12 else 1)
        assert p==expected_pass[(f,k)] and w==0
        rows.append(dict(family=f,support_weight=k,supports=s,
                         variants=v,linear_survivors=p,short_vectors=w))
    report=dict(
        claim="No vectors of norm <=8 with even first-half integer coordinates u, combined with existing layer4.py certificate",
        unproved="Odd-u case; full 96-dimensional minimum norm",
        A_sha256=candidate.get("A_sha256"),
        matrix_json_sha256=hashlib.sha256(json.dumps(A,separators=(',',':')).encode()).hexdigest(),
        binary_counts=dict(image_M={str(k):C0[k] for k in (11,12,14,16)},
                           image_Mt={str(k):C1[k] for k in (11,12,14,16)},
                           coset_ones_plus_M={str(k):C2[k] for k in (1,12,14,16)}),
        codewords_enumerated=3*(1<<24),finite_checks=rows,
        completed=True)
    Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
