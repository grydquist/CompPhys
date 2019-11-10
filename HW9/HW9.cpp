#include <iostream>
#include <functional>
#include "nr3.h"
#include "eigen_sym.h"

using namespace std;

//==============Funct(ions)/(tors)===============

// rhs
struct d {
    Int q;
    d(Int qq) : q(qq){};
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        dydx[0]= y[1];
        dydx[1]=(2*q*cos(2*x)-y[2])*y[0];
        dydx[2] = 0;
    }
};


int main(){
//======================P1=====================
Int N2 = 1;

for (Int i=0;i<15;i++){
    for (Int j=0;j<15;j++){
    };
}

//******************************************************************************

Int k = 20;
Doub q = 5.;

MatDoub a(k,k);

for (Int i=0; i<k; i++){
  for (Int j=0; j<k; j++){
    a[i][j]=0.;
    if (i==j){
      a[i][j] = 4.*i*j;
      if (j != k-1){
        a[i][j+1] = q;
      }
      if (j != 0){
        a[i][j-1] = q;
      }
    }

    }
  }


Jacobi  jac(a);

for (Int i=0; i<k; i++){
  cout <<a[i][i+1] << endl;
  cout << jac.d[i] << endl;
}

};
