// Exact rank-96 Construction-A norm <=8 search on first-half support six.
// Set oddmask bits for parity of the six nonzero coefficients (in support order).
// For each support and each of 4^6 compatible mod-8 residue assignments,
// use a necessary coordinatewise norm bound; enumerate EVERY centered
// mod-16 lift of each survivor; finish with exact integer squared norm.
// The first nonzero u coefficient is positive by global sign symmetry.
#include <immintrin.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix48.h"
using u64=unsigned long long;
static int oddmask=0, vodd=0;
static constexpr int sq16[16]={0,1,4,9,16,25,36,49,64,49,36,25,16,9,4,1};
static constexpr int sq8[8]={0,1,4,9,16,9,4,1};
static constexpr int odd[4]={1,3,5,7}, even[4]={0,2,4,6};
struct alignas(32) Vec {__m256i a,b;};
inline Vec plus(Vec x,Vec y){return {_mm256_add_epi8(x.a,y.a),_mm256_add_epi8(x.b,y.b)};}
inline Vec fromrow(int row,int coeff){
 alignas(32) uint8_t bytes[64]={};
 for(int j=0;j<48;j++)bytes[j]=(uint8_t)(coeff*A[row][j]);
 return {_mm256_load_si256((const __m256i*)bytes),
         _mm256_load_si256((const __m256i*)(bytes+32))};
}
inline int cost32(__m256i x,__m256i lut,__m256i mask){
 const __m256i s=_mm256_sad_epu8(
   _mm256_shuffle_epi8(lut,_mm256_and_si256(x,mask)),
   _mm256_setzero_si256());
 return (int)(_mm256_extract_epi64(s,0)+_mm256_extract_epi64(s,1)
             +_mm256_extract_epi64(s,2)+_mm256_extract_epi64(s,3));
}
struct Stats {
 u64 supports=0,classes=0,survive=0,lifts=0,shorts=0;
 int best=INT32_MAX,witness[48]={};
};
inline int residue(int pos,int digit){
 return ((oddmask>>pos)&1 ? odd : even)[digit];
}
inline int cost_u(int pos,int coeff){
 if(coeff==0)return 64; // support is nonzero: residue zero means +8 or -8
 if(pos==0)return coeff*coeff; // first nonzero coordinate has positive sign
 return sq8[coeff];
}
inline void check(int ids[6],Stats&z,__m256i lut,__m256i mask7){
 ++z.supports;
 Vec terms[6][4];
 for(int i=0;i<6;i++)for(int d=0;d<4;d++)
   terms[i][d]=fromrow(ids[i],residue(i,d));
 Vec left[64],right[64];
 int costleft[64],costright[64];
 for(int p=0;p<64;p++){
  int a=p>>4,b=(p>>2)&3,c=p&3;
  left[p]=plus(plus(terms[0][a],terms[1][b]),terms[2][c]);
  costleft[p]=cost_u(0,residue(0,a))+cost_u(1,residue(1,b))+cost_u(2,residue(2,c));
 }
 for(int q=0;q<64;q++){
  int a=q>>4,b=(q>>2)&3,c=q&3;
  right[q]=plus(plus(terms[3][a],terms[4][b]),terms[5][c]);
  costright[q]=cost_u(3,residue(3,a))+cost_u(4,residue(4,b))+cost_u(5,residue(5,c));
 }
 for(int p=0;p<64;p++)for(int q=0;q<64;q++){
  ++z.classes;
  int score=costleft[p]+costright[q];
  if(score+vodd>128)continue;
  Vec v=plus(left[p],right[q]);
  score+=cost32(v.a,lut,mask7);
  if(score>128)continue;
  score+=cost32(v.b,lut,mask7);
  if(score>128)continue;
  ++z.survive;
  int digit[6]={p>>4,(p>>2)&3,p&3,q>>4,(q>>2)&3,q&3};
  int r[6];for(int i=0;i<6;i++)r[i]=residue(i,digit[i]);
  int u[6];u[0]=(r[0]==0?8:r[0]);
  for(int bits=0;bits<32;bits++){
   for(int i=1;i<6;i++)
     u[i]=(r[i]?r[i]:8)-((bits>>(i-1)&1)?(r[i]?8:16):0);
   int sq=0;for(int i=0;i<6;i++)sq+=u[i]*u[i];
   if(sq+vodd>128)continue;
   ++z.lifts;
   for(int j=0;j<48 && sq<=128;j++){
    int w=0;for(int i=0;i<6;i++)w+=u[i]*A[ids[i]][j];
    sq+=sq16[((unsigned)w)&15];
   }
   if(sq<=128){
    ++z.shorts;
    if(sq<z.best){
     z.best=sq;std::memset(z.witness,0,sizeof(z.witness));
     for(int i=0;i<6;i++)z.witness[ids[i]]=u[i];
     #pragma omp critical
     {
      std::printf("WITNESS num=%d mask=%d support=%d,%d,%d,%d,%d,%d\n",
        sq,oddmask,ids[0],ids[1],ids[2],ids[3],ids[4],ids[5]);
      std::fflush(stdout);
     }
    }
   }
  }
 }
}
void recursive(int depth,int start,int ids[6],Stats&z,
               __m256i lut,__m256i mask7){
 if(depth==6){check(ids,z,lut,mask7);return;}
 for(int j=start;j<=48-(6-depth);j++){
  ids[depth]=j;
  recursive(depth+1,j+1,ids,z,lut,mask7);
 }
}
int main(int argc,char**argv){
 oddmask=argc>1?std::atoi(argv[1]):-1;
 const int begin=argc>2?std::atoi(argv[2]):0;
 const int end=argc>3?std::atoi(argv[3]):43;
 if(oddmask<0||oddmask>63||begin<0||end>43||end<=begin)return 2;
 const int k=__builtin_popcount((unsigned)oddmask);
 vodd=(k&1)?48-k:k; // A mod 2 = J-I
 alignas(32) uint8_t bytes[32];
 for(int i=0;i<32;i++)bytes[i]=sq8[i&7];
 const __m256i lut=_mm256_load_si256((const __m256i*)bytes),
               mask7=_mm256_set1_epi8(7);
 Stats total;
 #pragma omp parallel
 {
  Stats z;
  #pragma omp for schedule(dynamic,1)
  for(int first=begin;first<end;first++){
   int ids[6]={first};
   recursive(1,first+1,ids,z,lut,mask7);
  }
  #pragma omp critical
  {
   total.supports+=z.supports;total.classes+=z.classes;
   total.survive+=z.survive;total.lifts+=z.lifts;total.shorts+=z.shorts;
   if(z.best<total.best){
    total.best=z.best;
    std::memcpy(total.witness,z.witness,sizeof(z.witness));
   }
  }
 }
 std::printf(
  "SUMMARY oddmask=%d first=[%d,%d) supports=%llu classes=%llu survivors=%llu lifts=%llu short=%llu best_num=%d\n",
  oddmask,begin,end,total.supports,total.classes,total.survive,
  total.lifts,total.shorts,total.best);
 if(total.shorts){
  std::printf("WITNESS_U");
  for(int i=0;i<48;i++)std::printf(" %d",total.witness[i]);
  std::printf("\n");
 }
 return total.shorts?3:0;
}
