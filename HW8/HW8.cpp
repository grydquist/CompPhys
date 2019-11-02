#include <iostream>
#include <functional>
#include "nr3.h"
#include "stepper.h"
#include "odeint.h"
#include "stepperdopr853.h"
#include "shoot.h"
#include "shootf.h"
#include "ludcmp.h"
#include "qrdcmp.h"
#include "roots_multidim.h"

using namespace std;

//==============Funct(ions)/(tors)===============



Doub func(Doub x, Doub y){
    return exp(-x*sin(y));
};

struct d {
    Int q;
    d(Int qq) : q(qq){}; 
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        dydx[0]= y[1];
        dydx[1]=(2*q*cos(2*x)-y[2])*y[0];
        dydx[2] = 0;
    }
};

struct Load {
    Int n,q;
    VecDoub y;
    Load(Int nn, Int qq) : q(qq), n(nn), y(3) {};
    VecDoub operator() (const Doub x1, VecDoub_I &v) {
        y[2] = v[0];
        y[1] = 0;
        y[0] = 1;
        return y;
    }
};

struct Score {
    Int n,q;
    VecDoub f;
    Score(Int nn, Int qq) : n(nn), q(qq), f(1){};
    VecDoub operator() (const Doub xf, VecDoub_I &y){
        f[0] = y[0];
        return f;
    }
};

int main(){
//======================P1=====================
Int N2 = 1;
Int nvar = 3;
Int n = 3;
Int q=0;
VecDoub v(N2);
Bool check;

v[0] = n*n;

Load load(n,q);
d rhs(q);
Score score(n,q);

Doub x1 = 0.0;
Doub x2 = 3.1415926535897932384626433/2;
Shoot<Load,d,Score> shoot(nvar,x1,x2,load,rhs,score);
newt(v,check,shoot);
if (check) {
cout << "shoot failed; bad initial guess" << endl;
} else {
    cout<< v[0]<<endl;
};



};