#include <iostream>
#include <functional>
#include "nr3.h"
#include "eigen_sym.h"

using namespace std;

//==============Funct(ions)/(tors)===============

// rhs
struct rhs_func {
    Int q;
    rhs_func(Int qq) : q(qq){};
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        dydx[0]= y[1];
        dydx[1]=(2*q*cos(2*x)-y[2])*y[0];
    }
};


int main(){
//******************************************************************************
  cout.precision(12);

  Int k = 20;
  Doub q = 25.;

  MatDoub aodd(k,k), aeven(k,k), bodd(k,k), beven(k,k);

  for (Int i=0; i<k; i++){
    for (Int j=0; j<k; j++){
      aodd[i][j]=0.;
      beven[i][j]=0.;
    }
  }

  for (Int i=0; i<k; i++){
    for (Int j=0; j<k; j++){
      if (i==j){
        // aodd[i][j] = 4.*i*j;
        beven[i][j]=2.*(i+1)*2.*(j+1);
        aodd[i][j] = (2.*i+1.)*(2.*j+1.);

        if (j != k-1){
          aodd[i][j+1] = q;
          beven[i][j+1]=q;
        }
        if (j != 0){
          aodd[i][j-1] = q;
          beven[i][j-1]=q;
        }
      }
    }
  }
  bodd=aodd;
  aodd[0][0] = aodd[0][0] +q;
  bodd[0][0] = bodd[0][0] -q;


  Jacobi  jac_aodd(aodd), jac_bodd(bodd), jac_beven(beven);

  // eigsrt(jac.d);

  // for (Int i=0; i<k; i++){
  //   cout <<jac_aodd.d[i]<< endl;
  // }
  //
  // for (Int i=0; i<k; i++){
  //   cout <<jac_bodd.d[i]<< endl;
  // }

  cout<< endl
    <<"-------------------------------------------" << endl
    << left
    << setw(10) <<"n"
    <<"a,"<<" q="<<q<< endl
    <<"-------------------------------------------"<<endl
    << setw(10) <<"1"
    << setw(15)<<jac_aodd.d[k-1]<< endl
    << setw(10) <<"15"
    << setw(15)<< jac_aodd.d[k-8]<<endl;




  cout<< endl
    <<"-------------------------------------------" << endl
    << left
    << setw(10) <<"n"
    <<"b,"<<" q="<<q<< endl
    <<"-------------------------------------------"<<endl
    << setw(10) <<"1"
    << setw(15)<<jac_bodd.d[k-1]<< endl
    << setw(10) <<"2"
    << setw(15)<<jac_beven.d[k-1]<< endl
    << setw(10) <<"10"
    << setw(15)<< jac_beven.d[k-5]<<endl
    << setw(10) <<"15"
    << setw(15)<< jac_bodd.d[k-8]<<endl;

    //******************************************************************************
    Doub a = 7.449;
  return 0;
};
