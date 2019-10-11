#include <iostream>
#include <functional>
#include "nr3.h"
#include "ran.h"
#include "mcintegrate.h"

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

int main(){
//======================P2=====================
// box will be a square with bounds +- sqrt(5/8), MAYBE???
// achieved by minimizing ellipse eqn wrt x & y
Doub xlo,xhi,ylo,yhi;
xlo = -3;//-sqrt(5./8.);
ylo =  xlo;
xhi =  3;//sqrt(5./8.);
yhi =  xhi;

//MCintegrate myMC(xlo,xhi,vol,inell,1,1);

Doub n = 4e6;
Int xi = sqrt(n),yi = sqrt(n);
Doub xstep = (xhi-xlo)/xi,tot = 0,xv=xlo,yv=ylo;
Int xc=0,yc=0;
Bool in=false;

for (Int i=0;i<n;i++){
   // cout<<xv<<"    "<<yv<<endl;
    in = inell(xv,yv);
    if(in){
        tot++;
    }
    xv +=xstep;
    if (xv>xhi+xstep){
        xv = xlo;
        yv += xstep;
    }
};
cout<<tot<<endl;

tot *=(xhi-xlo)*(yhi-ylo)/n;

cout<<tot<<endl;

};