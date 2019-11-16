#include <iostream>
#include <functional>
#include <fstream>
#include "nr3.h"
#include "fourier.h"

using namespace std;

int main(){

  ofstream myfile;
  myfile.open("output.txt");
  myfile.precision(12);

  Doub a=1., V=1.;
  Doub pi = 3.1415926535897932384626433;
  Int nx=16, ny=nx;
  std::vector<Doub> an(nx);
  std::vector<std::vector<Doub> > u(nx,std::vector<Doub>(ny)),rho(nx,std::vector<Doub>(ny)),
  rhot(nx,std::vector<Doub>(ny)),uh(nx,std::vector<Doub>(ny));
  Doub del = 1/nx;

  for (Int i=0; i<nx; i++){
      an[i] = 1.;
    }

  sinft(an);
  
   for (Int i=0; i<nx; i++){
     an[i] *= 2./nx;
      // cout << an[i] << endl;
     }

  for (Int i=1; i<nx-1; i++){
      an[i] = 1./(sinh(pi*(i)))*an[i];
      //cout << an[i] << endl;
    }
    for (Int i=0; i<nx; i++){
      for (Int j=0; j<nx; j++){
        u[i][j]=0.;
        for (Int n=1; n<nx-1; n++){
          u[i][j]+= 2./nx * an[n]*sinh(pi*(n)*(i+1)/nx)*sin(pi*(n)*(j+1)/nx);
          }
        //cout << u[i][j] << endl;
        }
      }


//=================New stuff================
  std::vector<Doub> aa(nx);


// Make rho, where it's just 0 except -1 on the boundaries
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      rho[i][j]=0.;
      if(j == 1) rho[i][j] = -1;
    }
  }

// sine transform of rows (?)
  for (Int i=0; i<ny;i++){
      sinft(rho[i]);
  }

// swap indices so we can do sine transform on columns
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      rhot[i][j] = rho[j][i];
    }
  }
  rho = rhot;

// Sine transform on columns
  for (Int i=0; i<ny;i++){
      sinft(rho[i]);
  }

// swap rows back to original configuration
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      rhot[i][j] = rho[j][i];
    }
  }
  rho = rhot;

// solve for uh
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      uh[i][j] = rho[i][j]/2./(cos(pi*i/nx)+cos(pi*j/ny)-2);
    }
  }

// sine transform uh rows
  for (Int i=0; i<ny;i++){
      sinft(uh[i]);
  }

// swap rows and make above the inverse sine transform
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      u[i][j] = uh[j][i]*2/nx;
    }
  }
  uh = u;

// sine transform uh columns
  for (Int i=0; i<ny;i++){
      sinft(uh[i]);
  }

// Swap rows and columns again. Make above inverse
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      u[i][j] = uh[j][i]*2/nx;
    }
  }

// output the answer
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      if (abs(u[i][j])<1e-15) u[i][j]=0;
      std::cout << std::setprecision(5)<< u[i][j]<< " ";
      myfile    << std::setprecision(5)<< u[i][j]<< " ";
    }
      cout<<endl;
      myfile <<endl;
  }
   
  myfile.close(); 

    

  return 0;
};
