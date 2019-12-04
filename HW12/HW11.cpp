#include <iostream>
#include <functional>
#include "nr3.h"
#include "ludcmp.h"


using namespace std;

//==============Funct(ions)/(tors)===============

// Analytic solution
Doub fa(Doub x, Doub t) {
  return exp(-M_PI*M_PI*t)*sin(M_PI*x);
}

//=============Main Program====================
int main(){
  ofstream myfile;
  cout.precision(12);
// Declarations
  Doub L=2.,dt = 0.0001,t=0, t_end=1.;
  Int nx=40, counter=0;
  VecDoub x(nx),u(nx),ua(nx),ut(nx),ui(nx);
  Doub dx = L/(nx-1);
  Doub res_ftcs=0., res_cn=0.;

  x[0] = -1.;

// Setting x, analytic, numerical, temporary holder solutions
  for (Int i=0;i<nx;i++){
    if (i!=0){
    x[i] = x[i-1]+dx;
    }
    u[i] = fa(x[i],0.);
  }
  ua = u;
  ui = u;

// -----------CN--------------------
  MatDoub A(nx,nx);
  VecDoub b(nx);
  VecDoub t_array(round(1/dt));
  myfile.open("CN.txt");
  myfile.precision(12);

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

  A[0][0] = 1.;
  A[0][1] = 0.;
  A[1][0] = 0.;

  A[nx-1][nx-1] = 1.;
  A[nx-2][nx-1] = 0.;
  A[nx-1][nx-2] = 0.;

//Do LUdcmp of A
  LUdcmp ALU(A);

counter = 0;
//Time loop
  while (t<t_end){
    counter++;
//  Update time
    t = t+dt;
    t_array[counter]=t;
//  Loop through points to get RHS vector
    for (Int j=0;j<nx;j++){


//    Deal with points on ends of truncated domain
      if (j==0) {
        b[j] =0.;

      } else if (j==nx-1) {
        b[j] = 0.;

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
      myfile << std::setprecision(12)<< u[j]<<" ";
      res_cn+= (u[j]-ua[j])*(u[j]-ua[j]);
    }
    myfile << endl;
  }

  cout << sqrt(res_cn/counter/nx)<< endl;
  cout << counter << endl;
  myfile.close();

  myfile.open("xcn.txt");
  for (int i=0; i<nx; i++){
    myfile<<x[i]<<endl;
  }
  myfile.close();

  myfile.open("tcn.txt");
  for (int i=0; i<counter; i++){
    myfile<<t_array[i]<<endl;
  }
  return 0;
};
