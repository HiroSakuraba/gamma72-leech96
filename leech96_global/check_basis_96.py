#!/usr/bin/env python3
"""Verify original Z/16 Construction-A Gram and integral reduced basis exactly.

The determinant of the original 96x96 Gram is one by its Schur complement;
this script additionally checks the exported basis transform is unimodular
and its conjugated Gram equals the supplied reduced matrix.
"""
from pathlib import Path
import sys, json, hashlib, time
from sympy import Matrix

def matmul(X,Y):
    n=len(X); mid=len(Y); cols=len(Y[0]);out=[[0]*cols for _ in range(n)]
    for i in range(n):
        for k,a in enumerate(X[i]):
            if not a:continue
            for j,b in enumerate(Y[k]):out[i][j]+=a*b
    return out

def gram(A):
    n=len(A)
    I=[[int(i==j) for j in range(n)] for i in range(n)]
    AAt=matmul(A,list(map(list,zip(*A))))
    assert all((AAt[i][j]+int(i==j))%16==0 for i in range(n) for j in range(n))
    Q=[[(AAt[i][j]+int(i==j))//16 for j in range(n)] for i in range(n)]
    B=[[*(Q[i]),*(A[i])] for i in range(n)]
    B.extend([[*(A[k][i] for k in range(n)),*(16*int(i==j) for j in range(n))]
              for i in range(n)])
    assert all(B[i][j]==B[j][i] for i in range(2*n) for j in range(2*n))
    assert all(B[i][i]%2==0 for i in range(2*n))
    return B

def verify(candidate_path, transform_path, reduced_path):
    candidate=json.loads(Path(candidate_path).read_text())
    A=candidate['A'];n=len(A)
    assert n==48 and all(len(r)==n for r in A)
    h=hashlib.sha256(b''.join(int(x).to_bytes(8,'little',signed=True) for row in A for x in row)).hexdigest()
    assert h==candidate['A_sha256'],('candidate fingerprint mismatch',h)
    G=gram(A)
    U=json.loads(Path(transform_path).read_text())
    reduced=json.loads(Path(reduced_path).read_text())
    assert len(U)==2*n and all(len(r)==2*n for r in U)
    assert len(reduced)==2*n and all(len(r)==2*n for r in reduced)
    t=time.perf_counter()
    print('Computing unimodular transform determinant exactly...',flush=True)
    det=int(Matrix(U).det(method='domain-ge'))
    print('det(U)=',det,'seconds',round(time.perf_counter()-t,2),flush=True)
    assert abs(det)==1,'LLL transform fails unimodularity'
    T=matmul(matmul(U,G),list(map(list,zip(*U))))
    assert T==reduced,'reduced Gram differs from exact conjugated Gram'
    from exact_svp import exact_ldlt,verify_factor,digest
    piv,off=exact_ldlt(T)
    assert piv[-1]==1
    verify_factor(T,piv,off)
    report={
       'candidate_A_sha256':h,
       'reduced_gram_sha256':digest(reduced),
       'basis_transform_det':det,
       'even':all(reduced[i][i]%2==0 for i in range(2*n)),
       'reduced_gram_determinant':piv[-1],
       'reduced_GS_last_norm_squared':str(__import__('fractions').Fraction(piv[-1],piv[-2])),
       'basis_is_exactly_full_candidate_lattice':True,
       'shortest_vector_certified':False,
    }
    print(json.dumps(report,indent=2))
    return report

def main():
    if len(sys.argv) not in (4,5):raise SystemExit('Usage: check_basis_96.py candidate.json U.json G_reduced.json [report.json]')
    result=verify(*sys.argv[1:4])
    if len(sys.argv)==5:Path(sys.argv[4]).write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
