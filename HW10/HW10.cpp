#include <iostream>
#include <functional>
#include <fstream>
#include "nr3.h"
#include "fourier.h"

using namespace std;

//==============Funct(ions)/(tors)===============

// rhs
struct rhs_func {
    Int counter;
    Doub q, a;
    rhs_func(Doub aa, Doub qq) : a(aa), q(qq), counter(0){};
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        dydx[0]= y[1];
        dydx[1]=(2*q*cos(2*x)-a)*y[0];
    }
};


int main(){
  Doub a=1., V=1.;
  Doub pi = 3.1415926535897932384626433;
  Int nx=16, ny=nx;
  std::vector<Doub> an(nx);
  std::vector<std::vector<Doub> > u(nx,std::vector<Doub>(ny)),rho(nx,std::vector<Doub>(ny));
  Doub del = 1/nx;

  for (Int i=0; i<nx; i++){
      an[i] = 1.;
    }

  sinft(an);
  
   for (Int i=0; i<nx; i++){
     an[i] *= 2./nx;
      // cout << an[i] << endl;
     }

  for (Int i=1; i<nx-1; i++){
      an[i] = 1./(sinh(pi*(i)))*an[i];
      //cout << an[i] << endl;
    }
    for (Int i=0; i<nx; i++){
      for (Int j=0; j<nx; j++){
        u[i][j]=0.;
        for (Int n=1; n<nx-1; n++){
          u[i][j]+= 2./nx * an[n]*sinh(pi*(n)*(i+1)/nx)*sin(pi*(n)*(j+1)/nx);
          }
        //cout << u[i][j] << endl;
        }
      }


//=================New stuff================
    for (Int i=0; i<nx; i++){
      for (Int j=0; j<nx; j++){
        rho[i][j]=0.;
        if(j==nx-1) rho[i][j] = -1;
        }
      }

      for (Int i=0; i<ny;i++){
        sinft(rho[i]);
      }

    for (Int i=0; i<nx; i++){
      for (Int j=0; j<nx; j++){
        //cout<< rho[i][j]<<endl;
        u[i][j]=rho[i][j]/(cos(pi*i/nx)+cos(pi*j/ny)-2);
        }
      }

      for (Int i=0; i<ny;i++){
        sinft(u[i]);
      }

    for (Int i=0; i<nx; i++){
      for (Int j=0; j<nx; j++){
        //u[i][j] *= 2/nx;
        cout<< u[i][j]<<endl;
        }
      }

    

    

  return 0;
};
