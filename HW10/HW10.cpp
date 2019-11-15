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
  Int nx=16, ny=nx;
  std::vector<Doub> an(nx);
  std::vector<std::vector<Doub> > u(nx,std::vector<Doub>(ny));

  for (Int i=0; i<nx; i++){
      an[i] = 1.;
    }
  sinft(an);

  // sinft(an);
  // for (Int i=0; i<nx; i++){
  //     cout << an[i] * 2./nx << endl;
  //   }

  for (Int i=1; i<nx-1; i++){
      an[i] = 1./(sinh(M_PI*(i)/nx))*an[i];
      // cout << an[i] << endl;
    }
    for (Int i=0; i<nx; i++){
      for (Int j=0; j<nx; j++){
        u[i][j]=0.;
        for (Int n=1; n<nx-1; n++){
          u[i][j]+= 2./nx * an[n]*sinh(M_PI*(n)*(i+1)/nx)*sin(M_PI*(n)*(j+1)/nx);
          }
        cout << u[i][j] << endl;
        }
      }


  return 0;
};
