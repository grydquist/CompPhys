#include <iostream>
#include <functional>
#include "nr3.h"
#include "ludcmp.h"


using namespace std;

//==============Funct(ions)/(tors)===============

// Analytic solution
Doub fa(Doub x, Doub t) {
  Doub a = 1.+4.*t;
  return 1./pow(a,0.5)*exp(-1.*x*x/a);
}

//=============Main Program====================
int main(){
  ofstream myfile;
  myfile.open("ftcs.txt");
  myfile.precision(12);
// Declarations
  Doub L=10.,dt = 0.001,t=0;
  Int nx=100, Nt=1000;
  VecDoub x(nx),u(nx),ua(nx),ut(nx),ui(nx);
  Doub dx = L/(nx-1);
  x[0] = 0;

// Setting x, analytic, numerical, temporary holder solutions 
  for (Int i=0;i<nx;i++){
    if (i!=0){
    x[i] = x[i-1]+dx;
    }
    u[i] = fa(x[i],0.);
  }
  ua = u;
  ui = u;

// -------------FTCS----------------
// Time loop
  for (Int i=0;i<Nt;i++){
//  Update time
    t = t+dt;
//  Loop over all nodes to explicitly update next time step
    for (Int j=0;j<nx;j++){

//    Deal with points on ends of truncated domain
      if (j==0) {
        ut[j] = u[j] + dt*(2.*u[j+1]-2.*u[j])/dx/dx;

      } else if (j==nx-1) {
        ut[j] = u[j] + dt*(u[j-1]-2.*u[j])/dx/dx;

//    Update next timestep for interior points  
      } else {
        ut[j] = u[j] + dt*(u[j-1]-2.*u[j]+u[j+1])/dx/dx;
      }

//    Get analytical solution for current point/time
      ua[j] = fa(x[j],t);

//    Output info
      myfile << std::setprecision(5)<< ut[j]<<", "<<ua[j]<< endl;
    }

//  Update u to next time step
    u = ut;
    
  }
  myfile.close();
// -----------CN--------------------
  MatDoub A(nx,nx);
  VecDoub b(nx);

  myfile.open("CN.txt");
  t = 0;
  u = ui;
  ua = ui;

//Build LHS matrix, which is tridiagonal except the first off diagonal term
  for (Int i=0;i<nx;i++){
    for (Int j=0;j<nx;j++){
      A[i][j] = 0;
    }
  }

  for (Int i=0;i<nx;i++){
    A[i][i] = 1./dt + 1./dx/dx;
    if(i<nx-1) A[i+1][i] = -1./2./dx/dx;
    if(i>0)    A[i-1][i] = -1./2./dx/dx;
  }
  A[0][1] *=2;

//Do LUdcmp of A
  LUdcmp ALU(A);

//Time loop
  for (Int i=0;i<Nt;i++){

//  Update time
    t = t+dt;

//  Loop through points to get RHS vector
    for (Int j=0;j<nx;j++){


//    Deal with points on ends of truncated domain
      if (j==0) {
        b[j] = u[j]/dt+(2.*u[j+1]-2.*u[j])/2./dx/dx;

      } else if (j==nx-1) {
        b[j] = u[j]/dt+(u[j-1]-2.*u[j])/2./dx/dx;

//    Get RHS at points in interior
      } else {
        b[j] = u[j]/dt+(u[j-1]-2.*u[j]+u[j+1])/2./dx/dx;
      }

//    Analytical solution
      ua[j] = fa(x[j],t);
    }

//  Solve system
    ALU.solve(b,u);

//  Output
    for (Int j=0;j<nx;j++){
      myfile << std::setprecision(5)<< u[j]<<", "<<ua[j]<<", "<<x[j]<< endl;
    }
  }

  myfile.close();
  return 0;
};
