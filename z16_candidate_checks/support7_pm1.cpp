// Exhaustive exact test: u has 7 supported entries, each +/-1.
// First supported entry = +1 by global sign. 64 Gray-coded signs/support.
// This is an independent algorithmic cross-check of the larger all-odd scan.
#include <immintrin.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix48.h"
using u64=unsigned long long;
struct V {__m256i a,b;};
inline V add(V x,V y){return {_mm256_add_epi8(x.a,y.a),_mm256_add_epi8(x.b,y.b)};}
inline V sub(V x,V y){return {_mm256_sub_epi8(x.a,y.a),_mm256_sub_epi8(x.b,y.b)};}
inline V row(int i){return {_mm256_loadu_si256((const __m256i*)A_pad[i]),
                           _mm256_loadu_si256((const __m256i*)(A_pad[i]+32))};}
inline int cost(__m256i v,__m256i tab,__m256i mask){
 __m256i t=_mm256_shuffle_epi8(tab,_mm256_and_si256(v,mask));
 __m256i s=_mm256_sad_epu8(t,_mm256_setzero_si256());
 return (int)(_mm256_extract_epi64(s,0)+_mm256_extract_epi64(s,1)
             +_mm256_extract_epi64(s,2)+_mm256_extract_epi64(s,3));
}
struct Stats{
 u64 supports=0,signs=0,passed8=0,shorts=0;
 int best=INT_MAX,ids[7]={},sgn[7]={};
};
static __m256i tab8,tab16,m7,m15;
void check(const int*ids,Stats&st){
 ++st.supports;
 V rows[7];for(int i=0;i<7;i++)rows[i]=row(ids[i]);
 V v=rows[0];for(int i=1;i<7;i++)v=sub(v,rows[i]);
 int signs=0;
 for(int n=0;n<64;n++){
  ++st.signs;
  const int lower=7+cost(v.a,tab8,m7)+cost(v.b,tab8,m7);
  if(lower<=128){
   ++st.passed8;
   const int actual=7+cost(v.a,tab16,m15)+cost(v.b,tab16,m15);
   if(actual<=128){
    ++st.shorts;
    if(actual<st.best){
     st.best=actual;
     for(int j=0;j<7;j++){
      st.ids[j]=ids[j];
      st.sgn[j]=(j==0||((signs>>(j-1))&1))?1:-1;
     }
     #pragma omp critical
     {printf("SHORT numerator=%d first=%d\n",actual,ids[0]);fflush(stdout);}
    }
   }
  }
  if(n==63)break;
  const int bit=__builtin_ctz((unsigned)(n+1));
  const V twice=add(rows[bit+1],rows[bit+1]);
  if(!(signs&(1<<bit))){v=add(v,twice);signs|=(1<<bit);}
  else{v=sub(v,twice);signs&=~(1<<bit);}
 }
}
void rec(int at,int start,int*ids,Stats&st){
 if(at==7){check(ids,st);return;}
 for(int j=start;j<=48-(7-at);j++){
  ids[at]=j;rec(at+1,j+1,ids,st);
 }
}
int main(int argc,char**argv){
 int begin=argc>1?atoi(argv[1]):0,end=argc>2?atoi(argv[2]):42;
 if(begin<0||end>42||begin>=end)return 2;
 alignas(32)uint8_t a[32],b[32];
 const uint8_t s8[8]={0,1,4,9,16,9,4,1};
 const uint8_t s16[16]={0,1,4,9,16,25,36,49,64,49,36,25,16,9,4,1};
 for(int i=0;i<32;i++){a[i]=s8[i&7];b[i]=s16[i&15];}
 tab8=_mm256_load_si256((const __m256i*)a);
 tab16=_mm256_load_si256((const __m256i*)b);
 m7=_mm256_set1_epi8(7);m15=_mm256_set1_epi8(15);
 Stats total;
 #pragma omp parallel
 {
  Stats st;
  #pragma omp for schedule(dynamic,1)
  for(int i=begin;i<end;i++){int ids[7]={i};rec(1,i+1,ids,st);}
  #pragma omp critical
  {
   total.supports+=st.supports;total.signs+=st.signs;
   total.passed8+=st.passed8;total.shorts+=st.shorts;
   if(st.best<total.best){
    total.best=st.best;memcpy(total.ids,st.ids,sizeof(st.ids));
    memcpy(total.sgn,st.sgn,sizeof(st.sgn));
   }
  }
 }
 printf("SUMMARY first=[%d,%d) supports=%llu signed_patterns=%llu passed_mod8=%llu short=%llu best_numerator=%d\n",
        begin,end,total.supports,total.signs,total.passed8,total.shorts,total.best);
 if(total.shorts){
  printf("WITNESS");for(int j=0;j<7;j++)printf(" %d:%d",total.ids[j],total.sgn[j]);printf("\n");
 }
 return total.shorts?3:0;
}
