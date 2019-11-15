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
  Doub a=1.
  Int nx 

  return 0;
};
