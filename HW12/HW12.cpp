#include <iostream>
#include <functional>

#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
// #include "stepperbs.h"
#include "stepperdopr5.h"
#include "stepperdopr853.h"
#include "stepperbs.h"
#include "weights.h"
#include "ludcmp.h"
#include "stepperross.h"


using namespace std;

//==============Funct(ions)/(tors)===============

// Analytic solution
Doub fa(Doub x, Doub t) {
  return exp(-M_PI*M_PI*t)*sin(M_PI*x);
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

struct rhs_func2 {
  int counter, n;
  MatDoub D2;
  rhs_func2(MatDoub D22, Int nn) : counter(0), D2(D22), n(nn) {}
  void operator() (const Doub t, VecDoub_I &y, VecDoub_O &dydt){
    counter++;

    for (Int i=0; i<n/2+1; i++){
      dydt[i]=0.;
      for (Int j=0; j<n/2+1; j++){
        dydt[i]+=(D2[i][j]-D2[i][n-j-1])*y[j];
      }
    }
    dydt[0] = 0.;
    dydt[n/2]=0.;
  }
};

//=============Main Program====================
int main(){
  ofstream sol, analytical;
  sol.open("sol.txt");
  analytical.open("analytical.txt");
  sol.precision(12);
  analytical.precision(12);
  cout.precision(12);
// Declarations
  Doub atol=1.0e-15, rtol=atol, hw=0.01, hmin=0.0, t1=0.0, t2=1.0;
  Doub res=0.;
  Int n=21, nt=200, counter=0;
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
      // cout<< D2[i][j]<< " ";
    }
    // cout<< endl;
  }

  Output out_dopr(nt), out_dopr8(nt), out_bs(nt);
  rhs_func rhs(D2, n), rhs8(D2,n), rhsbs(D2,n);

  Odeint<StepperDopr5<rhs_func> > ode_dopr(ystart,t1,t2,atol,rtol,hw,hmin,out_dopr, rhs);
  Odeint<StepperDopr853<rhs_func> > ode_dopr8(ystart,t1,t2,atol,rtol,hw,hmin,out_dopr8, rhs8);
  // Odeint<StepperBS<rhs_func> > ode_bs(ystart,t1,t2,atol,rtol,hw,hmin,out_bs, rhsbs);

  ode_dopr.integrate();
  ode_dopr8.integrate();
  // ode_bs.integrate();
  // for (Int i=0; i<nt; i++){
  //   cout << out_dopr.ysave[1][i]<< " "<< fa(x[1],out_dopr.xsave[i])<<endl;
  // }

  VecDoub tsol(n);
  MatDoub ysol(n,nt+1);

  tsol = out_dopr.xsave;
  ysol = out_dopr.ysave;

  for (Int i=0; i<n; i++){
    for (Int j=0; j<=nt; j++){
      sol<< ysol[i][j]<<" ";
      analytical<< fa(x[i],tsol[j])<<" ";
      res+=(ysol[i][j]-fa(x[i],tsol[j]))*(ysol[i][j]-fa(x[i],tsol[j]));
    }
    sol << endl;
    analytical << endl;
  }
  cout<<"Number of rhs calls: "<<rhs.counter<<endl;

  cout<<sqrt(res/n/nt) << endl;

  sol.close();
  analytical.close();

  sol.open("time.txt");
  for (int j=0; j<=nt; j++){
    sol<< tsol[j]<< endl;
  }
  sol.close();

  sol.open("x_array.txt");
  for (int i=0; i<n; i++){
    sol<<x[i]<<endl;
  }
  sol.close();

  // Taking advantage of symmetry by only doing interval [0,1]
  VecDoub xi(n/2+1),ui(n/2+1);
  VecDoub ystarti(n/2+1);
  Doub res2=0.;
  MatDoub Ci(n/2+1,3), D2i(n/2+1,n/2+1);

  for (Int i=0;i<n/2+1 ;i++){
    xi[i] = sin(i*M_PI/(n-1));
    ystarti[i] = sin(M_PI*xi[i]);
    cout<<xi[i]<<endl;
  }
  for (Int i=0; i<n/2+1; i++){
    weights(xi[i],xi,Ci);
    for (Int j=0; j<n/2+1;j++){
      D2i[i][j]=Ci[j][2];
      // cout<< D2[i][j]<< " ";
    }
    // cout<< endl;
  }
  Output out_dopr_interval(nt);
  rhs_func rhs_interval(D2i, n/2+1);

  Odeint<StepperDopr5<rhs_func> > ode_dopr_interval(ystarti,t1,t2,atol,rtol,hw,hmin,out_dopr_interval, rhs_interval);

  ode_dopr_interval.integrate();

  VecDoub tsoli(n/2+1);
  MatDoub ysoli(n/2+1,nt+1);

  tsoli = out_dopr_interval.xsave;
  ysoli = out_dopr_interval.ysave;
  sol.open("interval.txt");

  for (Int i=0; i<n/2+1; i++){
    for (Int j=0; j<=nt; j++){
      sol<<ysoli[i][j]<<" ";
      res2+=(ysoli[i][j]-fa(xi[i],tsoli[j]))*(ysoli[i][j]-fa(xi[i],tsoli[j]));
    }
    sol<<endl;
  }
  cout<<sqrt(res2/nt/(n/2))<<endl;
  cout<<"Number of rhs calls: "<<rhs_interval.counter<<endl;
  sol.close();

  // Symmetry by only folding over D matrix, but keeping the interval [-1,1]
  VecDoub xs(n/2+1),us(n/2+1);
  VecDoub ystarts(n/2+1);
  Doub res3=0.;
  MatDoub Cs(n,3), D2s(n,n);


  for (Int i=0;i<n/2+1 ;i++){
    xs[i] = -cos(i*M_PI/(n-1));
    ystarts[i] = sin(M_PI*xs[i]);
    cout<<xs[i]<<" "<<ystarts[i]<<endl;
  }

  for (Int i=0; i<n; i++){
    weights(x[i],x,Cs);
    for (Int j=0; j<n;j++){
      D2s[i][j]=Cs[j][2];
    }
  }

  Output out_dopr_sym(nt);
  rhs_func2 rhs_sym(D2s, n);

  Odeint<StepperDopr5<rhs_func2> > ode_dopr_sym(ystarts,t1,t2,atol,rtol,hw,hmin,out_dopr_sym, rhs_sym);

  ode_dopr_sym.integrate();

  VecDoub tsols(n/2+1);
  MatDoub ysols(n/2+1,nt+1);

  tsols = out_dopr_sym.xsave;
  ysols = out_dopr_sym.ysave;
  sol.open("sym.txt");

  for (Int i=0; i<n/2+1; i++){
    for (Int j=0; j<=nt; j++){
      sol<<ysols[i][j]<<" ";
      res3+=(ysols[i][j]-fa(xs[i],tsols[j]))*(ysols[i][j]-fa(xs[i],tsols[j]));
    }
    sol<<endl;
  }

  cout << sqrt(res3/nt/(n/2+1))<<endl;

  return 0;
};
