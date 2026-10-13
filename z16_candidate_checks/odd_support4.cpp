// Exact support-four search for odd-u vectors in the published Z/16 lattice.
// Matrix A is generated and checked by odd_support4.py.
// u has exactly four nonzero centred mod-16 coordinates, v=uA (mod16).
#include <cstdint>
#include <cstdio>
#include <climits>
#include <cstdlib>
#include "matrix48.h"
#ifdef _OPENMP
#include <omp.h>
#endif
struct Stats {
  unsigned long long tested=0, parity=0, bound=0, complete=0;
  int best=INT_MAX;
  int u[48]={};
};
static constexpr int SQ[16]={0,1,4,9,16,25,36,49,64,49,36,25,16,9,4,1};
int main(int argc,char**argv){
  const int begin=argc>1?std::atoi(argv[1]):0;
  const int end=argc>2?std::atoi(argv[2]):45;
  const int maxabs=argc>3?std::atoi(argv[3]):8;
  if(begin<0 || end>45 || begin>=end || maxabs<1 || maxabs>8)return 2;
  Stats totals;
  #pragma omp parallel
  {
    Stats st;
    #pragma omp for schedule(dynamic,1)
    for(int i=begin;i<end;i++){
      for(int j=i+1;j<46;j++)for(int k=j+1;k<47;k++)for(int l=k+1;l<48;l++){
        for(int a=1;a<=maxabs;a++){
          int aa=a*a;
          int s0[48];for(int t=0;t<48;t++)s0[t]=a*A[i][t];
          for(int b=-maxabs;b<=maxabs;b++)if(b){
            int aabb=aa+b*b;
            if(aabb>=128)continue;
            int s1[48];for(int t=0;t<48;t++)s1[t]=s0[t]+b*A[j][t];
            for(int c=-maxabs;c<=maxabs;c++)if(c){
              int aabbcc=aabb+c*c;
              if(aabbcc>=128)continue;
              int s2[48];for(int t=0;t<48;t++)s2[t]=s1[t]+c*A[k][t];
              for(int d=-maxabs;d<=maxabs;d++)if(d){
                ++st.tested;
                int kk=(a&1)+(abs(b)&1)+(abs(c)&1)+(abs(d)&1);
                if(kk==0)continue;
                ++st.parity;
                int u2=aabbcc+d*d;
                if(u2+((kk&1)?48-kk:kk)>128)continue;
                ++st.bound;
                int cost=u2;
                for(int t=0;t<48;t++){
                  int residue=(s2[t]+d*A[l][t])&15;
                  cost+=SQ[residue];
                  if(cost>128)break;
                }
                if(cost<=128){
                  ++st.complete;
                  if(cost<st.best){
                    st.best=cost;st.u[i]=a;st.u[j]=b;st.u[k]=c;st.u[l]=d;
                    std::printf("FOUND numerator=%d i,j,k,l=%d,%d,%d,%d coeff=%d,%d,%d,%d\n",cost,i,j,k,l,a,b,c,d);
                    std::fflush(stdout);
                  }
                }
              }
            }
          }
        }
      }
    }
    #pragma omp critical
    {
      totals.tested+=st.tested;totals.parity+=st.parity;
      totals.bound+=st.bound;totals.complete+=st.complete;
      if(st.best<totals.best){totals.best=st.best;for(int i=0;i<48;i++)totals.u[i]=st.u[i];}
    }
  }
  std::printf("SUMMARY i=[%d,%d) maxabs=%d tested=%llu u_odd=%llu normbound=%llu short=%llu min_numerator=%d\n",begin,end,maxabs,totals.tested,totals.parity,totals.bound,totals.complete,totals.best);
  if(totals.complete){
    std::printf("WITNESS_U ");
    for(int i=0;i<48;i++)std::printf("%d%c",totals.u[i],i==47?'\n':' ');
  }
  return totals.complete?3:0;
}
