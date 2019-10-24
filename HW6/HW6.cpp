#include <iostream>
#include <functional>
#include "nr3.h"

using namespace std;

//==============Funct(ions)/(tors)===============

Bool inell(Doub x, Doub y){
    if (5.0*x*x-6.*x*y+5*y*y<2) {
        return true;
    }
    else{
        return false;
    }
};

Doub func(Doub x, Doub y){
    return exp(-x*sin(y));
}

int main(){
//======================P2=====================
Doub xlo,xhi,ylo,yhi;
xlo = -3;
ylo =  xlo;
xhi =  3;
yhi =  xhi;

Doub n = 9e4;
Int xi = sqrt(n),yi = sqrt(n), eval =0;
Doub xstep = (xhi-xlo)/xi,tot = 0,xv=xlo,yv=ylo;
Int xc=0,yc=0;
Bool in=false;

for (Int i=0;i<n;i++){
    in = inell(xv,yv);
    if(in){
        tot += func(xv,yv);
        eval++;
    }
    xv +=xstep;
    if (xv>xhi+xstep){
        xv = xlo;
        yv += xstep;
    }
};

tot *=(xhi-xlo)*(yhi-ylo)/n;
cout<< "Value of absolute error: " << abs(1.449026 - tot)<<endl;
cout<< "Number of integrand evaluations: " << eval<<endl;

};