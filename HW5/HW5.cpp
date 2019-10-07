#include <iostream>
#include <functional>
#include "nr3.h"
#include "dawson.h"
#include "chebyshev.h"

using namespace std;

//==============Funct(ions)/(tors)===============

struct func1{
    Doub E,El,Et;
    Doub x;

    func1(Doub xx) : x(xx) {};

    Doub next(Int n){
        Et = E;
        E = 1./n*(exp(-x)-x*El);
        El = Et;
        Et = 0.;
        if(n==1){
            El = 1;
            E = 0;
        }
    };
};

struct lentz{
    Doub eps=1;
    Doub D=0.,C, delt;
    Doub tiny = 1.e-30;
    Doub f;
    lentz(Doub f0) : C(f0),f(f0){
        if(f0==0) {C = tiny;f=tiny;};
    };
    Doub next(Doub a, Doub b){
        D = b + a*D;
        if(D==0) {D = tiny;};
        C = b+a/C;
        if(C==0) {C = tiny;};
        D = 1/D;
        delt = C*D;
        f = f*delt;
        eps = abs(delt-1);
    }
};

struct daws{
    Doub a,b,b0=0,x,n = 0;
    daws(Doub xx) : x(xx) {
        a = 0.5*xx;
        b = 0.5+xx*xx;
    };
    Doub operator() (Doub nn) {
        a = -(nn-1)*x*x;
        b = 0.5+x*x+(nn-1);
        n = n;
    };
};

Doub func2(Doub y) {
    return 1./(y*y*y-1) -1./(y*y*y-1.)/pow((y*y*y-1.),2./3.);
};

Doub func3 (Doub y) {
    return y*y*y-1;
}

Doub func4 (Doub x) {
    return 1./x - 1./x/pow(1+x,2./3.);
}





int main(){


//======================P1=====================
    func1 f1(1);

    //while (abs(f1.E-f1.El))
    for (Int n=1;n<20;n++) {
        f1.next(n);
    };


//======================P2=====================

    daws d1(0.1),d2(1.1);
    lentz l1(d1.b0),l2(d2.b0);
    Doub n=1;
    l1.next(d1.a,d1.b);
    l2.next(d2.a,d2.b);

    while (l1.eps>1.e-7) {
        n++;
        d1(n);
        l1.next(d1.a,d1.b);
    };
    cout<<endl<<"For x = 0.1, we get for Dawson's integral: "
    <<l1.f<< " in "<<n<<" iterations at a tolerance of 1e-7"<<endl;
    cout<<"Calculated with the dawson routine: "<<dawson(0.1)<<endl;

    n=1;
    while (l2.eps>1.e-7) {
        n++;
        d2(n);
        l2.next(d2.a,d2.b);
    };
    cout<<endl<<"For x = 0.1, we get for Dawson's integral: "
    <<l2.f<< " in "<<n<<" iterations at a tolerance of 1e-7"<<endl;
    cout<<"Calculated with the dawson routine: "<<dawson(1.1)<<endl;


//======================P3=====================

    Chebyshev mycheb(func2,pow(0.5,1/3),pow(3.,1/3), 100);

    cout <<endl<<func2(1)<<endl;



}