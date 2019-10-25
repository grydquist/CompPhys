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
    Doub eps=1.0e-16;
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        eps = x;
        if (x==0){eps=1.0e-16;};
        dydx[0]= y[1];
        dydx[1]=-1/eps*y[1]-(1-1/eps/eps)*y[0];
        c++;
    }
};

int main(){
//======================P1=====================

Doub j1 = Mybess.jn(1,1.0),j1p,j0 = Mybess.jn(0,1.0),atol1 = 1.0e-6,
h1 = 0.001,x1 =1.0,x2 = 30.0,rtol1 = atol1, j30 = Mybess.jn(1,30.0),
rtol2 = 1.0e-10, atol2 = rtol2;
j1p = j0 - j1;
VecDoub yst(2);
yst[0] = j1;
yst[1] = j1p;
Output out1(20),out2(20),out3(20),out4(20),out5(20);
rhs d;


Odeint<StepperDopr5<rhs> > myode(yst,x1,x2,atol1,rtol1,h1,0.0,out1,d);
myode.integrate();

cout <<endl<< "For RK5, 1e-6, number of function evals: "<<d.c<<endl;
cout << "Error of J1(30): "<<
abs(out1.ysave[0][20]-j30)/abs(j30)<<endl<<endl;

d.c = 0;
yst[0] = j1;
yst[1] = j1p;



Odeint<StepperBS<rhs> > myode2(yst,x1,x2,atol1,rtol1,h1,0.0,out2,d);
myode2.integrate();


cout << "For BS, 1e-6, number of function evals: "<<d.c<<endl;
cout << "Error of J1(30): "<<
abs(out2.ysave[0][20]-j30)/abs(j30)<<endl<<endl;

d.c = 0;
yst[0] = j1;
yst[1] = j1p;



Odeint<StepperDopr5<rhs> > myode3(yst,x1,x2,atol2,rtol2,h1,0.0,out3,d);
myode3.integrate();

cout << "For RK5, 1e-10, number of function evals: "<<d.c<<endl;
cout << "Error of J1(30): "<<
abs(out3.ysave[0][20]-j30)/abs(j30)<<endl<<endl;

d.c = 0;
yst[0] = j1;
yst[1] = j1p;


Odeint<StepperBS<rhs> > myode4(yst,x1,x2,atol2,rtol2,h1,0.0,out4,d);
myode4.integrate();

cout << "For BS, 1e-10, number of function evals: "<<d.c<<endl;
cout << "Error of J1(30): "<<
abs(out4.ysave[0][20]-j30)/abs(j30)<<endl<<endl;


//======================P2=====================
yst[0] = 0.0;
yst[1] = 0.5;
x1 = 0.0;
x2 = 1.0;
Odeint<StepperBS<rhs> > myode5(yst,x1,x2,atol1,rtol1,h1,0.0,out5,d);
myode5.integrate();

cout << "For BS, 1e-10, number of function evals: "<<d.c<<endl;
cout << "Error of J1(30): "<<
abs(out5.ysave[0][20]-j1)/abs(j1)<<endl<<endl;


};