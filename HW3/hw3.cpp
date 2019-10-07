#include <iostream>
#include "nr3.h"
#include "romberg.h"
#include "qgleg.h"
#include "quad2d.h"

using namespace std;

// All my functions to be integrated====================================

struct func1{
    //Counter
    Int c = 0;
    Doub operator() (Doub x){
        Doub I = 1/(1+x*x);
        c++;
        return I;
    };
};

struct func2{
    //Counter
    Int c = 0;
    Doub operator() (Doub x){
        Doub I = abs(x*x*x);
        c++;
        return I;
    };
};

struct func3{
    //Counter
    Int c = 0;
    Doub operator() (Doub x){
        Doub I = sqrt(x)/sin(x);
        c++;
        return I;
    };
};

Doub f4(Doub x) {
    Doub I = exp(-x)/sqrt(x);
    return I;
}
Doub f5(Doub x) {
    Doub I = sin(x)/(x*x);
    return I;
}
Doub f6(Doub x, Doub y) {
    Doub I=exp(-sqrt(2.)*(x-y)/2.*sin(sqrt(2.)/2*(x+y)));
    return I;
}

// Functions for the bounds on problem 3
Doub fb_2(Doub x){
    Doub I = sqrt(1-x*x)/2.;
    return I;
}
Doub fb_1(Doub x){
    Doub I = -sqrt(1-x*x)/2.;
    return I;
}

// Main Program=========================================================
int main()
{
Doub a=-1, b=1, gleps=1, tmp,tmp2;
Int N,N3=0,N6=0,N9=0;

// Problem 1============================================================

// Romberg integration
// Function 1
cout <<endl<< "Problem 1a function 1:"<<endl;
func1 f1;
func2 f2;
func3 f3;

qromb(f1,a,b,1.0e-3);
cout <<"Number of iterations to reach tolerance of "
<< 1e-03<<": " << f1.c<<endl;
f1.c = 0;

qromb(f1,a,b,1.0e-6);
cout <<"Number of iterations to reach tolerance of "
<< 1e-06<<": " << f1.c<<endl;
f1.c = 0;

qromb(f1,a,b,1.0e-9);
cout <<"Number of iterations to reach tolerance of "
<< 1e-09<<": " << f1.c<<endl;
f1.c = 0;

// Function 2
cout <<endl<< "Problem 1a function 2:"<<endl;

tmp = qromb(f2,a,b,1.0e-3);
cout <<"Number of iterations to reach tolerance of "
<< 1e-03<<": " << f2.c<<endl;
f2.c = 0;

tmp = qromb(f2,a,b,1.0e-6);
cout <<"Number of iterations to reach tolerance of "
<< 1e-06<<": " << f2.c<<endl;
f2.c = 0;

tmp = qromb(f2,a,b,1.0e-9);
cout <<"Number of iterations to reach tolerance of "
<< 1e-09<<": " << f2.c<<endl;
f2.c = 0;


// Gauss-Legendre
// Function 1
cout <<endl<< "Problem 1b function 1:"<< endl;
N = 1;
tmp = qgleg(f1,a,b,N);
while (gleps>1e-9) {
    N++;
    tmp2 = qgleg(f1,a,b,N);
    gleps = abs(tmp-tmp2)/tmp;
    tmp = tmp2;
    if (gleps<1e-3 && N3 == 0){
        cout<<"Number of points to acieve relative error of "
        <<1e-3<<": "<<N<<endl;
        N3=1;
    }
    if (gleps<1e-6 && N6 == 0){
        cout<<"Number of points to acieve relative error of "
        <<1e-6<<": "<<N<<endl;
        N6=1;
    }
    if (gleps<1e-9 && N9 == 0){
        cout<<"Number of points to acieve relative error of "
        <<1e-9<<": "<<N<<endl;
        N9=1;
    }
}

N3=0;N6=0;N9=0;

// Function 2
cout <<endl<< "Problem 1b function 2:"<< endl;
N = 1;
gleps = 1;
tmp = qgleg(f2,a,b,N);
while (gleps>1e-9) {
    N++;
    tmp2 = qgleg(f2,a,b,N);
    gleps = abs(tmp-tmp2)/tmp;
    tmp = tmp2;
    if (gleps<1e-3 && N3 == 0){
        cout<<"Number of points to acieve relative error of "
        <<1e-3<<": "<<N<<endl;
        N3=1;
    }
    if (gleps<1e-6 && N6 == 0){
        cout<<"Number of points to acieve relative error of "
        <<1e-6<<": "<<N<<endl;
        N6=1;
    }
    if (gleps<1e-9 && N9 == 0){
        cout<<"Number of points to acieve relative error of "
        <<1e-9<<": "<<N<<endl;
        N9=1;
    }
}

// Problem 2============================================================

// Function 3
cout <<endl<< "Problem 2a function 3:"<< endl;
Midpnt<Doub(Doub)> q3(f3,0,M_PI/2);
qromo(q3,1.e-7);

// Function 4
cout <<endl<< "Problem 2a function 4:"<< endl;
// split into two integrals, 0->1, 1-> inf
Midpnt<Doub(Doub)> q41(f4,0.0,10.0);
tmp = qromo(q41,1.e-7);
Midinf<Doub(Doub)> q42(f4,10.0,1.e99);
tmp = tmp+qromo(q42,1.e-7);

// Function 5
cout <<endl<< "Problem 2a function 5:"<< endl;
Midinf<Doub(Doub)> q5(f5,M_PI/2,1.e99);
tmp = qromo(q5,1.e-7);

// Problem 3============================================================

tmp = quad2d(f6,-1.0,1.0,fb_1,fb_2);
cout <<endl<<"Value of integration over ellipse: "<< tmp<<endl;

return 0;
}