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

// rhs
struct d {
    Int q;
    d(Int qq) : q(qq){}; 
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        dydx[0]= y[1];
        dydx[1]=(2*q*cos(2*x)-y[2])*y[0];
        dydx[2] = 0;
    }
};

// ab = 0, finding a, ab = 1 finding b: a == even, b == odd
struct Load {
    Int n, ab;
    VecDoub y;
    Load(Int nn) :  n(nn), y(3) {};
    VecDoub operator() (const Doub x1, VecDoub_I &v) {
        
        y[2] = v[0];

        if (ab ==0) { //a (even)
            y[0] = 1;
            y[1] = 0;
        }
        else {       //b (odd)
            y[0] = 0;
            y[1] = 1;
        }

        return y;
    }
};

// n even: antiperiodic. n odd: periodic
struct Score {
    Int n, ab;
    VecDoub f;
    Score(Int nn) : n(nn), f(1){};
    VecDoub operator() (const Doub xf, VecDoub_I &y){
        if ( n % 2 != 0) { // antiperiodic
            if (ab ==0) { // a (even)
                f[0] = y[0];
            }
            else { // antiperiodic odd
                f[0] = y[1];
            }
        }
        else {              //periodic
            if (ab ==0) { // a (even)
                f[0] = y[1];
            }
            else { // periodic odd
                f[0] = y[0];

            }
        }

        return f;
    }
};

int main(){
//======================P1=====================
Int N2 = 1;
Int nvar = 3;
Int n = 1;
Int q=1;
Int ii=-20;
Int outpts = 150;
Int cnta=0,cntb = 0,cnt;
VecDoub v(N2), nvals(5),qs(2), as(15),bs(14),btm(14),atm(15);
Bool check, tst;

nvals[0]=0;nvals[1]=1;nvals[2]=2;nvals[3]=10;nvals[4]=15;

Load load(n);
d rhs(q);
Score score(n);

Doub x1 = 0.0,min;
Doub x2 = 3.1415926535897932384626433/2;
Doub atol1 = 1.0e-6,h1 = 0.001 ,x2i = 4*x2,rtol1 = atol1;
Output out1(outpts);
VecDoub yst(3);


Shoot<Load,d,Score> shoot(nvar,x1,x2,load,rhs,score);


rhs.q = 5;
load.ab = 1;
score.ab = 1;
ii -= (cntb+1)%2;
load.n = ii;
score.n = ii;
v[0] = ii;
load(x1,v);
score(x1,v); 

// b
while (cntb<15) {
    // set guess, keeping odd or even (i.e. we're searching for a certain v)
    ii += 2;
    ii -= (cntb+1)%2;
    v[0] = ii;

    // now we can calculate b(cntb)
    newt(v,check,shoot);

    // Did it work?
    if (check) {
        cout << "shoot failed; bad initial guess" << endl;
    }
    // If so, check if we get something different
    else {
    // First ones we just need to plug right in
    // Othrewise, plug right in to see if it is different from any bs val

        // test all elements in b for new element
        tst = true;
        for (Int i=0;i<14;i++){
            if (abs(bs[i] - v[0]) <1e-3 ) {
                // not a new value
                tst = false;

            }
        }

        if (v[0]>240) tst = false; // too big

        if (cntb == 0 || cntb == 1 || tst) {
            
            bs[cntb] = v[0];
            cntb++;
            load.n = cntb+1;
            score.n = cntb+1;
            v[0] = ii;
            
            // set BCs to even or odd
            load(x1,v);
            score(x1,v); 
        }
    }
};

// now just need to sort them
btm = bs;
min = btm[0];
cnt = 0;
for (Int i=0;i<15;i++){
    for (Int j=0;j<15;j++){
        if (bs[j]<min) {
            min = bs[j];
            cnt = j;
        }
    }
    btm[i] = bs[cnt];
    bs[cnt] = 1000;
    min = bs[cnt];
}
bs = btm;
bs[14] = btm[14];
for (Int i=0;i<15;i++){
    cout<<"b"<<i+1<<" = "<< bs[i]<<endl;
}


//a
load.ab = 0;
score.ab = 0;
ii = -20;
load.n = ii;
score.n = ii;
v[0] = ii;
load(x1,v);
score(x1,v);

while (cnta<16) {
    // set guess, keeping odd or even (i.e. we're searching for a certain v)
    ii += 2;
    ii -= cnta%2;
    v[0] = ii;

    //calc a(cnta)
    newt(v,check,shoot);

    //worked?
    if (check) {
    cout << "shoot failed; bad initial guess" << endl;
    }
    // First ones we just need to plug right in
    // Othrewise, plug right in to see if it is different from last even/odd val
    else {        
        // test all elements in b for new element
        tst = true;
        for (Int i=0;i<16;i++){
            if (abs(as[i] - v[0]) <1e-3 ) {
                // not a new value
                tst = false;

            }
        }

        if (v[0]>240) tst = false; // too big

        if (cnta == 0 || cnta == 1 || tst) {
            as[cnta] = v[0];
            cnta++;
            load.n = cnta;
            score.n = cnta;
            
            // set BCs to even or odd
            load(x1,v);
            score(x1,v); 
        }
    };

};

// now just need to sort them
atm = as;
min = atm[0];
cnt = 0;
for (Int i=0;i<16;i++){
    for (Int j=0;j<16;j++){
        if (as[j]<min) {
            min = as[j];
            cnt = j;
        }
    }
    atm[i] = as[cnt];
    as[cnt] = 1000;
    min = as[cnt];
}
as = atm;
as[15] = atm[15];
for (Int i=0;i<16;i++){
    cout<<"a"<<i<<" = "<< as[i]<<endl;
}

//===================================
// I've calculated all the a's and b's, now I'm just very inefficiently integrating

//=======================================
// n == 0, q= 5, a
v[0] = as[0];
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode1(yst,x1,x2i,atol1,rtol1,h1,0.0,out1,rhs);
myode1.integrate();
ofstream myfile ("an0q5.txt");
for (Int i=0;i<outpts+1;i++){
    myfile <<out1.xsave[i]<<"     "<< out1.ysave[0][i]<<endl;
}
myfile.close();



//=======================================
// n == 1, q = 5, b
score.n = nvals[1];                                                        //change all these 
load.n = nvals[1];
load.ab = 1;
score.ab = 1;                                                              //change all these 
v[0] = bs[0];                                                              // change to a or b and number
Output out2(outpts);                                                       //change #     
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];                                                        
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode2(yst,x1,x2i,atol1,rtol1,h1,0.0,out2,rhs);//change # ode and out
myode2.integrate();                                                       //change #  
ofstream myfile2 ("bn1q5.txt");                                            //change name and myfile 
for (Int i=0;i<outpts+1;i++){
    myfile2 <<out2.xsave[i]<<"     "<< out2.ysave[0][i]<<endl;             //change #  and myfile
}
myfile2.close();                                                            //change myfile


//=======================================
// n == 1, q = 5, a
score.n = nvals[1];                                                        //change all these 
load.n = nvals[1];
load.ab = 0;
score.ab = 0;                                                              //change all these 
v[0] = as[1];                                                              // change to a or b & #
Output out3(outpts);                                                       //change #     
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];                                                   
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode3(yst,x1,x2i,atol1,rtol1,h1,0.0,out3,rhs);//change # ode and out
myode3.integrate();                                                       //change #  
ofstream myfile3 ("an1q5.txt");                                            //change name and myfile 
for (Int i=0;i<outpts+1;i++){
    myfile3 <<out3.xsave[i]<<"     "<< out3.ysave[0][i]<<endl;             //change #  and myfile
}
myfile3.close();                                                            //change myfile

//=======================================
// n == 2, q = 5, b
score.n = nvals[2];                                                        //change all these 
load.n = nvals[2];
load.ab = 1;
score.ab = 1;                                                              //change all these 
v[0] = bs[1];                                                              // change to a or b & #
Output out4(outpts);                                                       //change #     
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];                                                   
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode4(yst,x1,x2i,atol1,rtol1,h1,0.0,out4,rhs);//change # ode and out
myode4.integrate();                                                       //change #  
ofstream myfile4 ("bn2q5.txt");                                            //change name and myfile 
for (Int i=0;i<outpts+1;i++){
    myfile4 <<out4.xsave[i]<<"     "<< out4.ysave[0][i]<<endl;             //change out  and myfile
}
myfile4.close();                                                            //change myfile


//=======================================
// n == 2, q = 5, a
score.n = nvals[2];                                                        //change all these 
load.n = nvals[2];
load.ab = 0;
score.ab = 0;                                                              //change all these 
v[0] = as[2];                                                              // change to a or b & #
Output out5(outpts);                                                       //change #     
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];                                                   
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode5(yst,x1,x2i,atol1,rtol1,h1,0.0,out5,rhs);//change # ode and out
myode5.integrate();                                                       //change #  
ofstream myfile5 ("an2q5.txt");                                            //change name and myfile 
for (Int i=0;i<outpts+1;i++){
    myfile5 <<out5.xsave[i]<<"     "<< out5.ysave[0][i]<<endl;             //change out  and myfile
}
myfile5.close();                                                            //change myfile


//=======================================
// n == 10, q = 5, b
score.n = nvals[3];                                                        //change all these 
load.n = nvals[3];
load.ab = 1;
score.ab = 1;                                                              //change all these 
v[0] = bs[10];                                                              // change to a or b & #
Output out6(outpts);                                                       //change #     
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];                                                   
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode6(yst,x1,x2i,atol1,rtol1,h1,0.0,out6,rhs);//change # ode and out
myode6.integrate();                                                       //change #  
ofstream myfile6 ("bn10q5.txt");                                            //change name and myfile 
for (Int i=0;i<outpts+1;i++){
    myfile6 <<out6.xsave[i]<<"     "<< out6.ysave[0][i]<<endl;             //change out  and myfile
}
myfile6.close();                                                            //change myfile


//=======================================
// n == 10, q = 5, a
score.n = nvals[3];                                                        //change all these 
load.n = nvals[3];
load.ab = 0;
score.ab = 0;                                                              //change all these 
v[0] = as[10];                                                              // change to a or b & #
Output out7(outpts);                                                       //change #     
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];                                                   
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode7(yst,x1,x2i,atol1,rtol1,h1,0.0,out7,rhs);//change # ode and out
myode7.integrate();                                                       //change #  
ofstream myfile7 ("an10q5.txt");                                            //change name and myfile 
for (Int i=0;i<outpts+1;i++){
    myfile7 <<out7.xsave[i]<<"     "<< out7.ysave[0][i]<<endl;             //change out  and myfile
}
myfile7.close();                                                            //change myfile


//=======================================
// n == 15, q = 5, a
score.n = nvals[4];                                                        //change all these 
load.n = nvals[4];
load.ab = 0;
score.ab = 0;                                                              //change all these 
v[0] = as[15];                                                              // change to a or b & #
Output out8(outpts);                                                       //change #     
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];                                                   
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode8(yst,x1,x2i,atol1,rtol1,h1,0.0,out8,rhs);//change # ode and out
myode8.integrate();                                                       //change #  
ofstream myfile8 ("an15q5.txt");                                            //change name and myfile 
for (Int i=0;i<outpts+1;i++){
    myfile8 <<out8.xsave[i]<<"     "<< out8.ysave[0][i]<<endl;             //change out  and myfile
}
myfile8.close();                                                            //change myfile




//=======================================
// n == 15, q = 5, b
score.n = nvals[4];                                                        //change all these 
load.n = nvals[4];
load.ab = 1;
score.ab = 1;                                                              //change all these 
v[0] = bs[14];                                                              // change to a or b & #
Output out9(outpts);                                                       //change #     
load(x1,v);
yst[0] = load.y[0];
yst[1] = load.y[1];
yst[2] = v[0];                                                   
x2i = 4*x2;
x1 = 0;
Odeint<StepperDopr853<d> > myode9(yst,x1,x2i,atol1,rtol1,h1,0.0,out9,rhs);//change # ode and out
myode9.integrate();                                                       //change #  
ofstream myfile9 ("bn15q5.txt");                                            //change name and myfile 
for (Int i=0;i<outpts+1;i++){
    myfile9 <<out9.xsave[i]<<"     "<< out9.ysave[0][i]<<endl;             //change out  and myfile
}
myfile9.close();                                                            //change myfile




};