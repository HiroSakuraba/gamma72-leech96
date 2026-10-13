// Exact sparse-support radius-8 check for the published Z/16 rank-96 lattice.
// Compile -DSUPPORT=5 for all five-coordinate vectors, or -DSUPPORT=6 for all-odd six-coordinate vectors.
// matrix48.h is generated from candidate_3cycle.json by support_mod8_verify.py.
// All mathematical comparisons use integers; AVX2 only accelerates cost accumulation.
#include <immintrin.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix48.h"
#ifndef SUPPORT
#define SUPPORT 5
#endif
static_assert(SUPPORT==5 || SUPPORT==6,"support must be 5 or 6");
using u64=unsigned long long;
constexpr int M=SUPPORT==5?8:4;
constexpr int SPLIT=SUPPORT==5?2:3;
constexpr int TP=64, TT=SUPPORT==5?512:64;
static constexpr int sq16[16]={0,1,4,9,16,25,36,49,64,49,36,25,16,9,4,1};
static constexpr int sq8[8]={0,1,4,9,16,9,4,1};
static constexpr int residues8[8]={0,1,2,3,4,5,6,7};
static constexpr int residues_odd[4]={1,3,5,7};
struct alignas(32) Vec {__m256i a,b;};
inline Vec add8(const Vec&u,const Vec&v){return {_mm256_add_epi8(u.a,v.a),_mm256_add_epi8(u.b,v.b)};}
inline Vec fromrow(int row,int coeff){
 alignas(32) uint8_t z[64]={};
 for(int j=0;j<48;j++)z[j]=(uint8_t)(coeff*A[row][j]);
 return {_mm256_load_si256((const __m256i*)z),_mm256_load_si256((const __m256i*)(z+32))};
}
inline int costvec(__m256i x,__m256i lut,__m256i mask){
 __m256i r=_mm256_shuffle_epi8(lut,_mm256_and_si256(x,mask));
 __m256i s=_mm256_sad_epu8(r,_mm256_setzero_si256());
 return int(_mm256_extract_epi64(s,0)+_mm256_extract_epi64(s,1)
    +_mm256_extract_epi64(s,2)+_mm256_extract_epi64(s,3));
}
struct Stats {
 u64 supports=0,classes=0,passes=0,lifts=0,shorts=0;
 u64 passes_by_odd[7]={},lifts_by_odd[7]={};
 int best=999999,witness[48]={};
};
inline int residue(int i){return SUPPORT==5?residues8[i]:residues_odd[i];}
inline int min_u(int coeff,bool first){
 if(SUPPORT==5 && coeff==0)return 64; // u=+/-8: nonzero, even though residue=0
 return first?coeff*coeff:sq8[coeff];
}
void check_support(const int*idx,Stats &st,__m256i lut,__m256i mask7){
 ++st.supports;
 Vec terms[SUPPORT][M];
 for(int i=0;i<SUPPORT;i++)for(int j=0;j<M;j++)terms[i][j]=fromrow(idx[i],residue(j));
 Vec left[TP],right[TT];
 int lcost[TP],rcost[TT],lc[TP][SPLIT],rc[TT][SUPPORT-SPLIT];
 for(int p=0;p<TP;p++){
  int x=p,sum=0;
  Vec v={_mm256_setzero_si256(),_mm256_setzero_si256()};
  for(int a=SPLIT-1;a>=0;a--){
   int digit=x%M;x/=M;lc[p][a]=digit;
   v=add8(v,terms[a][digit]);
   sum+=min_u(residue(digit),a==0);
  }
  left[p]=v;lcost[p]=sum;
 }
 for(int q=0;q<TT;q++){
  int x=q,sum=0;
  Vec v={_mm256_setzero_si256(),_mm256_setzero_si256()};
  for(int a=SUPPORT-SPLIT-1;a>=0;a--){
   int digit=x%M;x/=M;rc[q][a]=digit;
   v=add8(v,terms[SPLIT+a][digit]);
   sum+=min_u(residue(digit),false);
  }
  right[q]=v;rcost[q]=sum;
 }
 for(int p=0;p<TP;p++)for(int q=0;q<TT;q++){
  ++st.classes;
  int lower=lcost[p]+rcost[q];
  if(lower>128)continue;
  Vec v=add8(left[p],right[q]);
  lower+=costvec(v.a,lut,mask7);
  if(lower>128)continue;
  lower+=costvec(v.b,lut,mask7);
  if(lower>128)continue;
  ++st.passes;
  int coeffs[SUPPORT];
  for(int a=0;a<SPLIT;a++)coeffs[a]=residue(lc[p][a]);
  for(int a=0;a<SUPPORT-SPLIT;a++)coeffs[SPLIT+a]=residue(rc[q][a]);
  int odd=0;for(int a=0;a<SUPPORT;a++)odd+=(coeffs[a]&1);
  ++st.passes_by_odd[odd];
  int u[SUPPORT];u[0]=(coeffs[0]?coeffs[0]:8);
  for(int lifts=0;lifts<(1<<(SUPPORT-1));lifts++){
    for(int a=1;a<SUPPORT;a++){
      int r=coeffs[a];u[a]=(r?r:8)-(((lifts>>(a-1))&1)?(r?8:16):0);
    }
    int numerator=0;for(int a=0;a<SUPPORT;a++)numerator+=u[a]*u[a];
    if(numerator+((odd&1)?48-odd:odd)>128)continue;
    ++st.lifts;++st.lifts_by_odd[odd];
    for(int col=0;col<48&&numerator<=128;col++){
      int w=0;for(int a=0;a<SUPPORT;a++)w+=u[a]*A[idx[a]][col];
      numerator+=sq16[((unsigned)w)&15];
    }
    if(numerator<=128){
      ++st.shorts;
      if(numerator<st.best){
        st.best=numerator;
        std::memset(st.witness,0,sizeof(st.witness));
        for(int a=0;a<SUPPORT;a++)st.witness[idx[a]]=u[a];
        #pragma omp critical
        {
         std::printf("WITNESS numerator=%d support=",numerator);
         for(int a=0;a<SUPPORT;a++)std::printf("%d%c",idx[a],a==SUPPORT-1?'\n':',');
         std::fflush(stdout);
        }
      }
    }
  }
 }
}
void enumerate(int depth,int start,int *ids,Stats&z,__m256i lut,__m256i mask7){
 if(depth==SUPPORT){check_support(ids,z,lut,mask7);return;}
 for(int j=start;j<=48-(SUPPORT-depth);j++){
   ids[depth]=j;
   enumerate(depth+1,j+1,ids,z,lut,mask7);
 }
}
int main(int argc,char**argv){
 const int begin=argc>1?std::atoi(argv[1]):0;
 const int end=argc>2?std::atoi(argv[2]):49-SUPPORT;
 if(begin<0||end>49-SUPPORT||end<=begin)return 2;
 alignas(32) uint8_t bytes[32];for(int i=0;i<32;i++)bytes[i]=(uint8_t)sq8[i&7];
 const __m256i lut=_mm256_load_si256((const __m256i*)bytes);
 const __m256i mask7=_mm256_set1_epi8(7);
 Stats total;
 #pragma omp parallel
 {
  Stats z;
  #pragma omp for schedule(dynamic,1)
  for(int first=begin;first<end;first++){
    int ids[SUPPORT]={};ids[0]=first;
    enumerate(1,first+1,ids,z,lut,mask7);
  }
  #pragma omp critical
  {
   total.supports+=z.supports;total.classes+=z.classes;
   total.passes+=z.passes;total.lifts+=z.lifts;total.shorts+=z.shorts;
   for(int j=0;j<=SUPPORT;j++){
    total.passes_by_odd[j]+=z.passes_by_odd[j];
    total.lifts_by_odd[j]+=z.lifts_by_odd[j];
   }
   if(z.best<total.best){total.best=z.best;std::memcpy(total.witness,z.witness,sizeof(z.witness));}
  }
 }
 std::printf("SUMMARY support=%d first=[%d,%d) supports=%llu mod8_classes=%llu pass_mod8=%llu lifts=%llu shorts=%llu best_num=%d\n",
 SUPPORT,begin,end,total.supports,total.classes,total.passes,total.lifts,total.shorts,total.best);
 std::printf("MOD8_PASSES_BY_ODD_COUNTS");for(int i=0;i<=SUPPORT;i++)std::printf(" %llu",total.passes_by_odd[i]);std::printf("\n");
 std::printf("LIFTS_BY_ODD_COUNTS");for(int i=0;i<=SUPPORT;i++)std::printf(" %llu",total.lifts_by_odd[i]);std::printf("\n");
 if(total.shorts){std::printf("WITNESS_U ");for(int i=0;i<48;i++)std::printf("%d%c",total.witness[i],i==47?'\n':' ');}
 return total.shorts?3:0;
}
