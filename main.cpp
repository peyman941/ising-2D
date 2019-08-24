#include "iostream"
#include "fstream"
#include "ctime"
#include "cstdlib"
#include "math.h"
///////////////////////////////////
using namespace std;
double ran2(long *idum);
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
long SEED = -1973;
//
// Paul de Bakker, 6 Aug 2000
//
// FROM NUMERICAL RECIPES IN C
//
// Long period (>2 E18) random number generator of
// L'Ecuyer with Bays-Durham suffle and added
// safeguards. Returns a uniform random deviate between
// 0.0 and 1.0 (exclusive of the endpoint values).
// Call with idum a negative integer to initialise;
// thereafter, do not alter idum between successive
// deviates in a sequence. RNMX should approximate the
// largest floating value that is less than 1.
//
double ran2(long *idum)
{
int j;
long k;
static long idum2 = 123456789;
static long iy = 0;
static long iv[NTAB];
double temp;
if (*idum <= 0)
{
if (-(*idum) < 1)
*idum = 1;
else *idum = -(*idum);
idum2 = (*idum);
for (j = NTAB+7; j >= 0; j--)
{
k = (*idum) / IQ1;
*idum = IA1 * (*idum - k*IQ1) - k*IR1;
if (*idum < 0)
*idum += IM1;
if (j < NTAB)
iv[j] = *idum;
}
iy = iv[0];
}
k = (*idum) / IQ1;
*idum = IA1 * (*idum - k*IQ1) - k*IR1;if (*idum < 0)
*idum += IM1;
k = idum2 / IQ2;
idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
if (idum2 < 0)
idum2 += IM2;
j = iy / NDIV;
iy = iv[j] - idum2;
iv[j] = *idum;
if (iy < 1)
iy += IMM1;
if ((temp = AM * iy) > RNMX)
return RNMX;
else
return temp;
} // ran2
int main() {
    ofstream output("./ising23.txt");
    int N=100;
    //int T0=0;
    double T=0.0;
    double* energyt=new double[60];
    double* meght=new double[60];
    double* energyt2=new double[60];
    double* meght2=new double[60];

    for(int i=0; i<60; i++) {
      energyt[i]=0.0;
      meght[i]=0.0;
      energyt2[i]=0.0;
      meght2[i]=0.0;
    }

    double avgenergy0,energy0,energy,avgenergy;//energy
    double megh0,avgmegh0,avgmegh,megh;//meghnatesh
    energy0=0;
    megh0=0;
    int net[N][N];
    //unsigned seed=time(NULL);
    //srand(seed);
    //initial condition we set +1 to all the spin


  /*  std::cout << net[i][j] << '\n';*/

/*for (size_t i = 0; i <100; i++) {
for (size_t j = 0; j <100; j++) {
     energy0+=-(net[i][j])*(net[(i-1)%N][j]+net[(i+1)%N][j]+net[i][(j+1)%N]+net[i][(j-1)%N]);
     megh0+=net[i][j];
  }
}

  avgenergy0=(energy0/10000);
  avgmegh0=(megh0/10000);
    output<<T0<<"\t\t"<<avgenergy0<<"\t\t"<<avgmegh0<<"\t\t"<<endl;
    energy=energy0;
    megh=megh0;*/
    //////////////////////////////////////////////
    for (size_t av = 0; av < 100; av++) {

      for (size_t i = 0; i <100; i++) {
      for (size_t j = 0; j <100; j++) {
      net[i][j]=(+1);
    }
  } // in each avarage we set the spin up(it oue tirick)

  for (int t= 0;t <60; t++) {
    long seed =32.0; //seed of my random number
  for (size_t n=0; n<1000000; n++) {
     T=0.1*t;
     int c=ran2(&seed)*100;
     int d=ran2(&seed)*100;
     energy=2*((net[c][d])*(net[(c-1)%N][d]+net[(c+1)%N][d]+net[c][(d+1)%N]+net[c][(d-1)%N]));
     if (energy<0) {
       net[c][d]=-net[c][d];
      }
      if(energy>0){
       double a=ran2(&seed);
        if(a<exp((-1.0)*energy/(T+0.0))){
          net[c][d]=-net[c][d];
         }
       }
   }
   //energyt[t]=0.0;
    //meght[t]=0.0;
double energ=0.0,megh=0.0;
   for (size_t i = 0; i <100; i++) {
   for (size_t j = 0; j <100; j++) {

energ +=(-net[i][j])*(net[(i-1)%N][j]+net[(i+1)%N][j]+net[i][(j+1)%N]+net[i][(j-1)%N])/2.0;//specific energy
 megh +=net[i][j]; //calculate meghmatesh
}
}
energyt[t] +=energ; //total energy
 meght[t] +=megh;//tatoal  meghnatesh
 energyt2[t] +=energ*energ; //power 2 on energy
  meght2[t] +=megh*megh; //power 2 on meghnatesh
}
}
//avgenergy=(energyt[t]/10000);
//avgmegh=(meght[t]/10000);
for(int i=0;i<60;i++){  //because we have 5 temprature step with 0.1 arrival so we have 50 step
  energyt[i]=energyt[i]/100.0;  //avg on energy for each temp step
  meght[i]=meght[i]/100.0;  //avg on megh for each temp
  energyt2[i]=energyt2[i]/100.0;    //avg on power 2 energy
  meght2[i]=meght2[i]/100.0;  //avg on power 2 megh
output<<i*0.1<<"\t\t"<<energyt[i]/10000.0<<"\t\t"<<meght[i]/10000.0<<"\t\t"<<(meght2[i]-meght[i]*meght[i])/(1000*i)<<"\t\t"<<(energyt2[i]-energyt[i]*energyt[i])/(100*i*i)<<endl;
}
 return 0;
   }
