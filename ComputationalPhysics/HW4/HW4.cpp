#include <iostream>
#include <functional>
#include "nr3.h"
#include "roots.h"
#include "roots_poly.h"
#include "interp_1d.h"
#include "series.h"

using namespace std;

//==============Funct(ions)/(tors)===============


struct func1{
    Doub p;
    Doub operator() (const Doub x){
        return cos(x)-0.8+p*x*x;
    };
    Doub df(const Doub x) {
        return -sin(x)+2*p*x;
    };
};



struct func2{
    Doub x;
    Doub operator() (Doub p) {
        func1 f1;
        f1.p = p;
        using std::placeholders::_1;
        std::function<Doub(Doub)> df1 = std::bind(&func1::df,f1,_1);
        x = rtbis<std::function<Doub(Doub)>> (df1,0.01,10,1e-9);
        return f1(x);
    }
};


//==============Main Program=====================
int main()
{
// ----------------P1---------------------------
    func1 f1;
    func2 f2;
    VecDoub pp(3);
    Doub proot;
    for (Int i=1;i<4;i++) {
        f1.p = 0.1*i;

//      Four roots, but the function is symmetric about the y-axis
        cout <<endl<< "Roots for p = "<<f1.p<<":"<<endl;
        cout << rtnewt(f1,-5,-2,1e-3)<<endl;
        cout << rtnewt(f1,-2,0,1e-3)<<endl;
        cout << rtnewt(f1,0,2,1e-3)<<endl;
        cout << rtnewt(f1,2,5,1e-3)<<endl;
        }

        proot = rtbis<func2> (f2,0,0.4,1e-9);
        cout <<endl<<"There is a double root when p = "<<proot<<endl;
        cout<<"This double root occurs at x = "<<f2.x<<endl;


// -------------P2-------------------------
    VecComplex a(7),b(6);
    a[0] = 1;
    a[1] = -12.1;
    a[2] = 59.5;
    a[3] = -151.85;
    a[4] = 212.6625;
    a[5] = -156.6;
    a[6] = 48.5625;

    zroots(a,b,false);
    cout <<endl<< "Roots of the polynomial: "<<endl;
    for(Int i = 0;i<6;i++){
        cout <<b[i]<<endl;
    }

// ---------------P3-----------------------
    Eulsum myeul(100,1.e-09);
}
