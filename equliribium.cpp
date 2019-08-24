/*ising without periodic condition*/
#include "iostream"
#include "fstream"
#include "ctime"
#include "cstdlib"
#include "cmath"
using namespace std;

int main(int argc, char const *argv[]) {
double megh;
double energy;
    ofstream output("./ising11.txt");
    int N=100;
    //int T;
  int net[N][N];


  unsigned seed=time(NULL);
  srand(seed);

  for (size_t i = 0; i <100; i++) {
    for (size_t j = 0; j <100; j++) {
    net[i][j]=(+1);
    }
  }

  float T;
  for (size_t t = 0; t< 60.0; t++) {
 for (size_t n = 0; n<1000; n++) {
   T=0.1*t;

    int c=rand()%10;
    int d=rand()%10;

    energy=(net[c][d])*(net[(c-1)%N][d]+net[(c+1)%N][d]+net[c][(d+1)%N]+net[c][(d-1)%N]);

   if (energy<0) {
     net[c][d]=(-1)*net[c][d];
   }
   else if(energy<0){
     int a=(random()%10)/10;
     
     if(a<exp((2*energy)/T)){
        net[c][d]=(-1)*net[c][d];
   }
   
  }


  for (size_t i = 0; i <100; i++) {
    for (size_t j = 0; j <100; j++) {
megh+=net[i][j];
megh=megh/10000;
  energy+=energy;
  energy=energy/10000;
    }
  }
}
output<<"\t"<<megh<<"\t"<<"\t"<<T<<endl;

}
return 0;
}
