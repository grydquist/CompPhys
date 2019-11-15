#include <iostream>
#include <functional>
#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "eigen_sym.h"
#include "stepperdopr5.h"
#include <fstream>

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
  ofstream myfile;
  myfile.open("output.txt");
  myfile.precision(12);
//******************************************************************************
  cout.precision(16);

  Int k = 15;
  Doub q = 5.;

  MatDoub aodd(k,k), aeven(k,k), bodd(k,k), beven(k,k);

  for (Int i=0; i<k; i++){
    for (Int j=0; j<k; j++){
      aodd[i][j]=0.;
      beven[i][j]=0.;
    }
  }
  aeven=aodd;
  for (Int i=0; i<k; i++){
    for (Int j=0; j<k; j++){
      if (i==j){
        // aodd[i][j] = 4.*i*j;
        beven[i][j]=2.*(i+1)*2.*(j+1);
        aodd[i][j] = (2.*i+1.)*(2.*j+1.);
        aeven[i][j] = (2.*i)*(2.*j);

        if (j != k-1){
          aodd[i][j+1] = q;
          beven[i][j+1]=q;
          aeven[i][j+1] = q;
        }
        if (j != 0){
          aodd[i][j-1] = q;
          beven[i][j-1]=q;
          aeven[i][j-1] = q;
        }
      }
    }
  }
  bodd=aodd;
  aodd[0][0] = aodd[0][0] +q;
  bodd[0][0] = bodd[0][0] -q;

  aeven[0][1]=sqrt(2)*q;
  aeven[1][0]=sqrt(2)*q;

  Jacobi  jac_aodd(aodd), jac_aeven(aeven), jac_bodd(bodd), jac_beven(beven);

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
    << setw(10) <<"0"
    << setw(15)<<jac_aeven.d[k-1]<< endl
    << setw(10) <<"1"
    << setw(15)<< jac_aodd.d[k-1]<<endl
    << setw(10) <<"2"
    << setw(15)<<jac_aeven.d[k-2]<< endl
    << setw(10) <<"10"
    << setw(15)<< jac_aeven.d[k-6]<<endl
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
    cout<< endl
      <<"-------------------------------------------" << endl
      << " Integration from 0 to 2Pi" << endl;
      <<"-------------------------------------------"<<endl;

    Doub a = jac_aeven.d[k-2];
    const Int nvar = 2;
    Doub atol=1.0e-10, rtol=atol, hw=0.01, hmin=0.0, x1=0., x2=2.*M_PI;
    VecDoub ystart(nvar);

    ystart[0] = 1.;
    ystart[1] = 0.;

    rhs_func d1(a, q);
    Output out_dopr(50);
    Odeint<StepperDopr5<rhs_func> > ode_dopr(ystart,x1,x2,atol,rtol,hw,hmin,out_dopr,d1);

    ode_dopr.integrate();


    for (Int i=0; i<out_dopr.count; i++){
      cout << out_dopr.xsave[i]<<" "<<out_dopr.ysave[0][i]<<endl;
      myfile << out_dopr.xsave[i]<<" "<<out_dopr.ysave[0][i]<<endl;
    }
    myfile.close();

  return 0;
};
