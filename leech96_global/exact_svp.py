#!/usr/bin/env python3
"""Exact Fincke-Pohst enumeration for integral positive-definite Gram matrices.

An entirely integer/rational verifier with restartable deterministic DFS. No
floating-point arithmetic is used to prune search nodes. Enumeration searches
for nonzero x in Z^n with x^T G x <= bound, modulo the sign symmetry x~-x.

A completed witness-free search for a proven unimodular integral Gram basis
certifies absence of vectors of norm <= bound. A node-limited run NEVER does.

Usage:
    python exact_svp.py selftest
    python exact_svp.py inspect reduced_gram.json
    python exact_svp.py prefixes reduced_gram.json --depth 1 --output suffixes.json
    python exact_svp.py search reduced_gram.json --bound 8 --max-nodes 10000 --state run.json
    python exact_svp.py search reduced_gram.json --bound 8 --suffix 0 --max-nodes 10000 --state run0.json
    python exact_svp.py search reduced_gram.json --resume run.json --max-nodes 10000

The checkpoint contains its own stack plus exact rational remaining norms;
resume recomputes and checks all active choices. Certification requires the
Gram input and independent basis check, plus full completion of all suffixes.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import sys

VERSION=1

def digest(matrix):
    return sha256(json.dumps(matrix, separators=(',',':')).encode()).hexdigest()

def read_gram(path):
    G=json.loads(Path(path).read_text())
    n=len(G)
    if n==0 or any(not isinstance(row,list) or len(row)!=n for row in G):
        raise ValueError('Gram must be a square array')
    for i in range(n):
        for j in range(n):
            if type(G[i][j]) is not int or G[i][j]!=G[j][i]:
                raise ValueError('Gram must be symmetric integral')
    return G

def exact_ldlt(G):
    """Fraction-free symmetric Bareiss factorization.

    For each pivot i, p_i=det(G[:i+1,:i+1]), p_-1=1.
    The lower-triangular off[j][i] and pivot p_i satisfy
        x G x^t = sum_{i=0}^{n-1}
          (p_i*x_i + sum_{j>i} off[j][i]*x_j)^2 / (p_{i-1}*p_i).
    Exact division checks prove the elimination identity.
    """
    H=[row.copy() for row in G]
    n=len(G)
    off=[[0]*n for _ in range(n)]
    pivots=[]
    prev=1
    for k in range(n):
        pivot=H[k][k]
        if pivot<=0:raise ValueError(f'not positive definite, pivot {k}={pivot}')
        pivots.append(pivot)
        for j in range(k+1,n):off[j][k]=H[j][k]
        for i in range(k+1,n):
            a=H[i][k]
            for j in range(i,n):
                num=pivot*H[i][j]-a*H[j][k]
                val,rem=divmod(num,prev)
                if rem:raise AssertionError('nonintegral Bareiss division')
                H[i][j]=H[j][i]=val
        prev=pivot
    return pivots,off

def verify_factor(G,p,off):
    """Independently reconstruct the Gram using rational arithmetic.
    The Cauchy sums are exact and detect corrupted LDL factor data.
    """
    n=len(G)
    for i in range(n):
        for j in range(i,n):
            z=Fraction(0)
            for k in range(i+1):
                ai=p[k] if i==k else off[i][k]
                aj=p[k] if j==k else off[j][k]
                prev=p[k-1] if k else 1
                z+=Fraction(ai*aj,prev*p[k])
            if z!=G[i][j]:raise AssertionError((i,j,z,G[i][j]))

def nearest_first(low,high,shift,pivot):
    # Signs and x=0 first when convenient, then increasing exact squared GS cost.
    return sorted(range(low,high+1),key=lambda x:(abs(pivot*x+shift),abs(x),x))

def ceildiv(a,b):return -((-a)//b)

class Enumerator:
    def __init__(self, G, bound=8, suffix=()):
        self.G=G
        self.n=len(G)
        self.bound=int(bound)
        if self.bound<0:raise ValueError('bound must be nonnegative')
        self.matrix_hash=digest(G)
        self.p,self.off=exact_ldlt(G)
        self.shift=[0]*self.n
        self.x=[0]*self.n
        self.frames=[]
        self.nodes_total=0
        self.leaves_total=0
        self.witness=None
        self.suffix=list(suffix)
        self.finished=False
        self.seen=False
        self.start_i=self.n-1
        self.start_R=Fraction(bound)
        self._assign_suffix(self.suffix)
        if self.start_i<0:
            self.finished=True
            if self.seen:
                self._check_witness()
        else:
            self._push(self.start_i,self.start_R,self.seen)
            if not self.frames:self.finished=True

    def _apply(self,i,val):
        if self.x[i]!=0:raise AssertionError('already assigned')
        self.x[i]=val
        if val:
            for k in range(i):
                self.shift[k]+=self.off[i][k]*val

    def _undo(self,i):
        val=self.x[i]
        if val:
            for k in range(i):self.shift[k]-=self.off[i][k]*val
        self.x[i]=0

    def _bounds(self,i,R,seen):
        if R<0:raise AssertionError('negative remaining bound')
        w=(self.p[i-1] if i else 1)*self.p[i]
        lim=isqrt((R.numerator*w)//R.denominator)
        lo=ceildiv(-lim-self.shift[i],self.p[i])
        hi=(lim-self.shift[i])//self.p[i]
        if not seen:lo=max(lo,0)
        return lo,hi

    def _term(self,i,val):
        prev=self.p[i-1] if i else 1
        numer=(self.p[i]*val+self.shift[i])**2
        return Fraction(numer,prev*self.p[i])

    def _assign_suffix(self,suffix):
        if len(suffix)>self.n:raise ValueError('suffix longer than rank')
        R=Fraction(self.bound)
        seen=False
        for depth,val in enumerate(suffix):
            i=self.n-1-depth
            if type(val) is not int:raise ValueError('suffix must contain integers')
            lo,hi=self._bounds(i,R,seen)
            if not lo<=val<=hi:raise ValueError(f'inadmissible suffix at level {i}: {val} not in [{lo},{hi}]')
            term=self._term(i,val)
            R-=term
            self._apply(i,val)
            seen=seen or bool(val)
        self.seen=seen
        self.start_i=self.n-len(suffix)-1
        self.start_R=R

    def _push(self,i,R,seen):
        lo,hi=self._bounds(i,R,seen)
        if lo>hi:return False
        order=nearest_first(lo,hi,self.shift[i],self.p[i])
        self.frames.append({'i':i,'lo':lo,'hi':hi,'next':0,
                            'rn':R.numerator,'rd':R.denominator,
                            'seen':seen,'active':None,'order':order})
        return True

    def _check_witness(self):
        y=sum(self.x[i]*sum(self.G[i][j]*self.x[j] for j in range(self.n))
              for i in range(self.n))
        if not 0<y<=self.bound:
            raise AssertionError(f'pruned radius discrepancy at leaf: norm {y}')
        self.witness={'coefficients':self.x.copy(),'norm':y}

    def run(self, max_nodes=10000, stop_on_witness=True):
        max_nodes=int(max_nodes)
        if max_nodes<0:raise ValueError('node budget must be nonnegative')
        processed=0
        while self.frames and processed<max_nodes and not (stop_on_witness and self.witness):
            f=self.frames[-1]
            i=f['i']
            if f['active'] is not None:
                self._undo(i)
                f['active']=None
            if f['next']>=len(f['order']):
                self.frames.pop()
                continue
            val=f['order'][f['next']]
            f['next']+=1
            processed+=1
            self.nodes_total+=1
            R=Fraction(f['rn'],f['rd'])
            rem=R-self._term(i,val)
            if rem<0:raise AssertionError('pruning error')
            self._apply(i,val)
            f['active']=val
            seen=f['seen'] or bool(val)
            if i==0:
                self.leaves_total+=1
                if seen:self._check_witness()
                continue
            self._push(i-1,rem,seen)
        if not self.frames:self.finished=True
        return {'processed':processed,'complete':self.finished,
                'witness':self.witness is not None,
                'nodes_total':self.nodes_total,'leaves_total':self.leaves_total}

    def save(self,path):
        payload={'format':'exact-svp-state','version':VERSION,'gram_sha256':self.matrix_hash,
                 'rank':self.n,'bound':self.bound,'suffix':self.suffix,
                 'nodes_total':self.nodes_total,'leaves_total':self.leaves_total,
                 'witness':self.witness,'finished':self.finished,
                 'stack':[{k:v for k,v in f.items() if k!='order'} for f in self.frames],
                 'search_scope':'full lattice up to sign, restricted to suffix' if self.suffix else 'full lattice up to sign',
                 'certificate_status':('WITNESS' if self.witness else ('COMPLETE_NO_SHORT_VECTOR' if self.finished else 'INCOMPLETE')),
                 'global_minimum_certified':bool(self.finished and not self.witness and not self.suffix)}
        Path(path).write_text(json.dumps(payload,indent=2)+'\n')
        return payload

    @classmethod
    def resume(cls,G,payload):
        if payload['format']!='exact-svp-state' or payload['version']!=VERSION:
            raise ValueError('unsupported checkpoint version')
        if payload['gram_sha256']!=digest(G):raise ValueError('Gram checksum changed')
        obj=cls(G,payload['bound'],tuple(payload['suffix']))
        obj.nodes_total=payload['nodes_total']
        obj.leaves_total=payload['leaves_total']
        obj.witness=payload['witness']
        obj.finished=payload['finished']
        if obj.witness is not None:
            w=obj.witness
            z=w.get('coefficients')
            if not isinstance(z,list) or len(z)!=obj.n or not all(type(v) is int for v in z):
                raise ValueError('invalid witness in checkpoint')
            norm=sum(z[i]*sum(G[i][j]*z[j] for j in range(obj.n)) for i in range(obj.n))
            if norm!=w.get('norm') or not (0<norm<=obj.bound):
                raise ValueError('incorrect witness in checkpoint')
            if obj.suffix and [z[obj.n-1-j] for j in range(len(obj.suffix))]!=obj.suffix:
                raise ValueError('witness has wrong suffix')
        # Replay completed frame choices, checking exact remaining norms, bounds
        obj.frames=[]
        R=obj.start_R
        seen=obj.seen
        for idx,raw in enumerate(payload['stack']):
            i=obj.start_i-idx
            if raw['i']!=i or (raw['rn'],raw['rd'])!=(R.numerator,R.denominator) or raw['seen']!=seen:
                raise ValueError('inconsistent checkpoint frame')
            lo,hi=obj._bounds(i,R,seen)
            order=nearest_first(lo,hi,obj.shift[i],obj.p[i])
            if (lo,hi)!=(raw['lo'],raw['hi']) or not 0<=raw['next']<=len(order):
                raise ValueError('checkpoint bounds/iterator mismatch')
            active=raw['active']
            if active is not None:
                if raw['next']==0 or order[raw['next']-1]!=active:
                    raise ValueError('invalid active candidate')
            f=dict(raw,order=order)
            obj.frames.append(f)
            if active is None:
                if idx+1<len(payload['stack']):raise ValueError('unassigned parent')
            else:
                R-=obj._term(i,active)
                obj._apply(i,active)
                seen=seen or bool(active)
        if obj.finished and obj.frames:
            raise ValueError('completed search with nonempty stack')
        if not obj.finished and not obj.frames and not obj.witness:
            raise ValueError('missing pending search stack')
        return obj


def prefixes(G,depth,bound=8):
    """Exhaustively list sign-normalized admissible top-coordinate suffixes."""
    obj=Enumerator(G,bound)
    if depth<0 or depth>len(G):raise ValueError('bad prefix depth')
    out=[]
    def go(i,R,seen,chosen):
        if len(chosen)==depth:
            out.append(chosen.copy())
            return
        lo,hi=obj._bounds(i,R,seen)
        for val in range(lo,hi+1):
            term=obj._term(i,val)
            obj._apply(i,val)
            go(i-1,R-term,seen or bool(val),chosen+[val])
            obj._undo(i)
    go(len(G)-1,Fraction(bound),False,[])
    return out

def selftest():
    # E8 Cartan matrix, determinant 1, 240 roots of norm two.
    n=8
    G=[[2 if i==j else 0 for j in range(n)] for i in range(n)]
    for i,j in [(k,k+1) for k in range(6)]+[(2,7)]:G[i][j]=G[j][i]=-1
    p,o=exact_ldlt(G)
    assert p[-1]==1
    verify_factor(G,p,o)
    root=Enumerator(G,bound=1)
    assert root.run(max_nodes=1000000)['complete']
    assert root.witness is None and root.finished
    # E8 has exactly 240 roots of norm 2, hence 120 representatives
    # under x ~ -x, plus the zero vector.
    roots=Enumerator(G,bound=2)
    assert roots.run(max_nodes=1000000, stop_on_witness=False)['complete']
    assert roots.leaves_total==121 and roots.witness['norm']==2,roots.leaves_total
    assert len(prefixes(G,1,bound=2))>0
    # Basis transform catches rational pivots not equal to 1.
    A=[[2,-1],[-1,2]]
    verify_factor(A,*exact_ldlt(A))
    e=Enumerator(A,bound=1)
    assert e.run(100000)['complete'] and e.witness is None
    # Test checkpoint/resume bitwise equality of exact DFS results.
    Q=[[2,1],[1,2]]
    e=Enumerator(Q,bound=1)
    r=e.run(1)
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        fn=Path(td)/'checkpoint.json'
        data=e.save(fn)
        resumed=Enumerator.resume(Q,json.loads(fn.read_text()))
        resumed.run(100000)
        baseline=Enumerator(Q,bound=1)
        baseline.run(100000)
        assert resumed.nodes_total==baseline.nodes_total
        assert resumed.witness==baseline.witness
        assert resumed.finished==baseline.finished
    print('PASS exact Bareiss LDL, E8 norm-1 exclusion, suffix partition, checkpoint/resume')

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser('selftest')
    insp=sub.add_parser('inspect');insp.add_argument('gram')
    pr=sub.add_parser('prefixes');pr.add_argument('gram');pr.add_argument('--depth',type=int,default=1);pr.add_argument('--bound',type=int,default=8);pr.add_argument('--output',required=True)
    run=sub.add_parser('search');run.add_argument('gram');run.add_argument('--bound',type=int,default=8);run.add_argument('--max-nodes',type=int,default=10000);run.add_argument('--suffix',default='');run.add_argument('--state',default='svp_state.json');run.add_argument('--resume',default='')
    a=p.parse_args()
    if a.cmd=='selftest':selftest();return
    G=read_gram(a.gram)
    if a.cmd=='inspect':
        piv,off=exact_ldlt(G)
        verify_factor(G,piv,off)
        print(json.dumps({'rank':len(G),'gram_sha256':digest(G),'determinant':piv[-1],
              'first_pivots':piv[:4],'last_pivots':piv[-4:],
              'last_GS_norm_sq':[str(Fraction(piv[i],piv[i-1] if i else 1)) for i in range(len(G)-8,len(G))],
              'bound8_top_branch_count':len(prefixes(G,1,8))},indent=2))
        return
    if a.cmd=='prefixes':
        suff=prefixes(G,a.depth,a.bound)
        out={'status':'COMPLETE_PREFIX_LIST','gram_sha256':digest(G),'rank':len(G),
             'bound':a.bound,'depth':a.depth,'suffixes':suff}
        Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
        print('COMPLETE prefix list depth',a.depth,'count',len(suff))
        return
    if a.cmd=='search':
        if a.resume:
            data=json.loads(Path(a.resume).read_text())
            obj=Enumerator.resume(G,data)
        else:
            suffix=tuple(int(v) for v in a.suffix.split(',') if v.strip())
            obj=Enumerator(G,a.bound,suffix=suffix)
        result=obj.run(a.max_nodes)
        payload=obj.save(a.state)
        print(json.dumps({'processed_in_call':result['processed'],'nodes_total':obj.nodes_total,
                          'leaves_total':obj.leaves_total,'remaining_frames':len(obj.frames),
                          'status':payload['certificate_status'],
                          'witness':obj.witness,'checkpoint':a.state},indent=2))
        return

if __name__=='__main__':main()
