#include <iostream>
#include <functional>
#include "nr3.h"
#include "bessel.h"
#include "stepper.h"
#include "stepperdopr5.h"
#include "stepperbs.h"
#include "odeint.h"

using namespace std;

//==============Funct(ions)/(tors)===============
Bessjy Mybess;



Doub func(Doub x, Doub y){
    return exp(-x*sin(y));
};

struct rhs {
    Int c=0;
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        dydx[0]= y[1];
        dydx[1]=-1/x*y[1]-(1-1/x/x)*y[0];
        c++;
    }
};

int main(){
//======================P1=====================

Doub j1 = Mybess.jn(1,1.0),j1p,j0 = Mybess.jn(0,1.0),atol = 1.0e-6,
h1 = 0.01,x1 =1.0,x2 = 30.0,rtol = atol, j30 = Mybess.jn(1,30.0);
j1p = j0 - j1;
VecDoub yst(2);
yst[0] = j1;
yst[1] = j1p;
Output out(20);
rhs d;

Odeint<StepperDopr5<rhs> > myode(yst,x1,x2,atol,rtol,h1,0.0,out,d);
myode.integrate();

cout << "For RK5, number of function evals: "<<d.c<<endl;
cout << "Error of J1(30): "<<
abs(out.ysave[0][20]-j30)/j30<<endl;

d.c = 0;

Odeint<StepperBS<rhs> > myode(yst,x1,x2,atol,rtol,h1,0.0,out,d);

};