// Exact finite search of even-u cases in the published Z/16 candidate.
// Header matrix48.h is generated and validated by even_layer_certificate.py.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include "matrix48.h"
using u64=uint64_t;
static constexpr u64 FULL=(1ULL<<48)-1;
struct Counts {unsigned long long supports=0, cases=0, modular=0, short_vectors=0;};

int basis(const u64 *rows,u64 *out) {
  u64 piv[48]={};
  for(int i=0;i<48;i++){
    u64 x=rows[i];
    while(x){int p=63-__builtin_clzll(x);
      if(piv[p])x^=piv[p];else{piv[p]=x;break;}
    }
  }
  int r=0;
  for(int i=0;i<48;i++)if(piv[i])out[r++]=piv[i];
  return r;
}

// Necessary conditions modulo 8, solved over F_2.
// The sum(sign bits) mod 4 is relaxed; this permits false positives only.
// extra=-1 means no modification; otherwise add 2*a to w[extra].
// For extra outside support, t[extra] must be 2 or 6 modulo 8.
bool modular_possible(u64 supp,int k,int extra,int a,int desired){
  int S[16],J[48],ns=0,nj=0;
  for(int i=0;i<48;i++)
    if((supp>>i)&1)S[ns++]=i;else J[nj++]=i;
  if(ns!=k || nj!=48-k)std::abort();
  bool inside=extra>=0&&((supp>>extra)&1);
  for(int h=0;h<4;h++){
    u64 piv[48]={};
    int rhs_piv[48]={};
    bool impossible=false;
    for(int jj=0;jj<nj;jj++){
      int j=J[jj],sum=0;
      for(int s=0;s<ns;s++)sum+=A[S[s]][j];
      if(extra>=0)sum+=2*a*A[extra][j];
      int goal=(extra>=0&&!inside&&j==extra)?desired:0;
      if((sum-goal)&1)std::abort();
      int r=(((sum-goal)/2)%4+4)%4;
      int d=(h-r+4)%4;
      if(d&1){impossible=true;break;}
      int rhs=d>>1;
      u64 eq=column_masks[j]&supp;
      while(eq){int p=63-__builtin_clzll(eq);
        if(piv[p]){eq^=piv[p];rhs^=rhs_piv[p];}
        else{piv[p]=eq;rhs_piv[p]=rhs;break;}
      }
      if(!eq&&rhs){impossible=true;break;}
    }
    if(impossible)continue;
    if(inside){
      int rhs=a==1?0:1;
      u64 eq=1ULL<<extra;
      while(eq){int p=63-__builtin_clzll(eq);
        if(piv[p]){eq^=piv[p];rhs^=rhs_piv[p];}
        else{piv[p]=eq;rhs_piv[p]=rhs;break;}
      }
      if(!eq&&rhs)continue;
    }
    return true;
  }
  return false;
}

// Exhaust ALL 2^k sign patterns after the necessary linear filter.
// Norm=(|w|^2+|t|^2)/4 with t=wA mod8; centered t gives minimal norm.
bool actual_short(u64 supp,int k,int extra,int a){
  int S[16],ns=0;
  for(int i=0;i<48;i++)if((supp>>i)&1)S[ns++]=i;
  if(ns!=k)std::abort();
  bool inside=extra>=0&&((supp>>extra)&1);
  int special=-1;
  if(inside)for(int p=0;p<k;p++)if(S[p]==extra)special=p;
  int t[48]={};
  for(int s=0;s<k;s++)
    for(int j=0;j<48;j++)t[j]+=A[S[s]][j];
  if(extra>=0)for(int j=0;j<48;j++)t[j]+=2*a*A[extra][j];
  const int sq[8]={0,1,4,9,16,9,4,1};
  for(int step=0;step<(1<<k);step++){
    int gray=step^(step>>1);
    if(!inside||((gray>>special)&1)==(a==-1)){
      int cost=k+(extra<0?0:(inside?8:4));
      for(int j=0;j<48&&cost<=32;j++)cost+=sq[(t[j]%8+8)%8];
      if(cost<=32){
        std::printf("SHORT support=%llu k=%d extra=%d a=%d gray=%d numerator=%d\n",
                    (unsigned long long)supp,k,extra,a,gray,cost);
        return true;
      }
    }
    if(step+1==(1<<k))break;
    int flip=__builtin_ctz((unsigned)(step+1));
    int diff=((gray>>flip)&1)?2:-2;
    for(int j=0;j<48;j++)t[j]+=diff*A[S[flip]][j];
  }
  return false;
}

int main(){
  for(int family=0;family<3;family++){
    u64 B[48]={};
    int r=basis(family==1?column_masks:row_masks,B);
    if(r!=24)std::abort();
    for(int k:{12,16}){
      Counts c;
      u64 x=family==2?FULL:0;
      for(u64 n=1;n<(1ULL<<r);n++){
        x^=B[__builtin_ctzll(n)];
        if(__builtin_popcountll(x)!=k)continue;
        ++c.supports;
        ++c.cases;
        if(modular_possible(x,k,-1,0,0)){
          ++c.modular;
          if(actual_short(x,k,-1,0))++c.short_vectors;
        }
        if(k!=12)continue;
        for(int extra=0;extra<48;extra++){
          bool inside=(x>>extra)&1;
          for(int a:{-1,1}){
            ++c.cases;
            bool possible=inside?modular_possible(x,k,extra,a,0)
                  :(modular_possible(x,k,extra,a,2)||
                    modular_possible(x,k,extra,a,6));
            if(!possible)continue;
            ++c.modular;
            if(actual_short(x,k,extra,a))++c.short_vectors;
          }
        }
      }
      std::printf("family=%d k=%d supports=%llu cases=%llu modular_pass=%llu short=%llu\n",
        family,k,c.supports,c.cases,c.modular,c.short_vectors);
      std::fflush(stdout);
      if(c.short_vectors)return 2;
    }
  }
}
