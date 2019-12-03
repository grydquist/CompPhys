#include <iostream>
#include <functional>

#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperbs.h"
#include "stepperdopr5.h"
#include "stepperdopr853.h"
#include "weights.h"

using namespace std;

//==============Funct(ions)/(tors)===============

// Analytic solution
Doub fa(Doub x, Doub t) {
  Doub a = 1.+4.*t;
  return 1./pow(a,0.5)*exp(-1.*x*x/a);
}

struct rhs_func {
  int counter, n;
  MatDoub D2;
  rhs_func(MatDoub D22, Int nn) : counter(0), D2(D22), n(nn) {}
  void operator() (const Doub t, VecDoub_I &y, VecDoub_O &dydt){
    counter++;

    for (Int i=0; i<n; i++){
      dydt[i]=0.;
      for (Int j=0; j<n; j++){
        dydt[i]+=D2[i][j]*y[j];
      }
    }
    dydt[0] = 0.;
    dydt[n-1]=0.;
  }
};


//=============Main Program====================
int main(){
  ofstream myfile;
  myfile.open("ftcs.txt");
  myfile.precision(12);
  cout.precision(12);
// Declarations
  Doub atol=1.0e-6, rtol=atol, hw=0.01, hmin=0.0, t1=0.0, t2=1.0;
  Int n=10, counter=0;
  VecDoub x(n),u(n);
  VecDoub ystart(n);
  MatDoub C(n,3), D2(n,n);

  for (Int i=0;i<n ;i++){
    x[i] = -cos(i*M_PI/(n-1));
    ystart[i] = sin(M_PI*x[i]);
    cout<<x[i]<<endl;
  }

  for (Int i=0; i<n; i++){
    weights(x[i],x,C);
    for (Int j=0; j<n;j++){
      D2[i][j]=C[j][2];
      cout<< D2[i][j]<< " ";
    }
    cout<< endl;
  }

  Output out_dopr(20);
  rhs_func rhs(D2, n);

  Odeint<StepperDopr5<rhs_func> > ode_dopr(ystart,t1,t2,atol,rtol,hw,hmin,out_dopr, rhs);

  ode_dopr.integrate();

  for (Int i=0; i<20; i++){
    cout << out_dopr.ysave[1][i]<< " "<< exp(-M_PI*M_PI*out_dopr.xsave[i])*sin(M_PI*x[1])<<endl;
  }

  myfile.close();
  return 0;
};
