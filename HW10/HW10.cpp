#include <iostream>
#include <functional>
#include <fstream>
#include "nr3.h"
#include "fourier.h"

using namespace std;

int main(){

  ofstream myfile;
  myfile.open("sin.txt");
  myfile.precision(12);
  cout.precision(4);

  Doub a=1., V=1.;
  Doub pi = 3.1415926535897932384626433;
  Int nx=16, ny=nx;
  std::vector<Doub> an(nx);
  std::vector<std::vector<Doub> > u(nx,std::vector<Doub>(ny)),rho(nx,std::vector<Doub>(ny)),
  rhot(nx,std::vector<Doub>(ny)),uh(nx,std::vector<Doub>(ny));
  Doub del = 1./nx;

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
  cout<< endl<<
  "-----------------------------------------------------------------------------" <<endl
  <<"Sine Transform" << endl<<
  "-----------------------------------------------------------------------------"
  <<endl;
// output the answer
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      if (abs(u[i][j])<1e-15) u[i][j]=0;
      cout<< left<< setw(8)<<u[i][j];
      myfile    << std::setprecision(5)<< u[i][j]<< " ";
    }
      cout<<endl;
      myfile <<endl;
  }

  myfile.close();

  //=================Cosine Tranform attempt ================
  myfile.open("cos.txt");
  nx=17; ny=nx;
  del = 1./nx;

  std::vector<std::vector<Doub> > uc(nx,std::vector<Doub>(ny)),rhoc(nx,std::vector<Doub>(ny)),
  rhotc(nx,std::vector<Doub>(ny)),uhc(nx,std::vector<Doub>(ny));

// Make rho, where it's just 0 except -1 on the boundaries
cout<<
"-----------------------------------------------------------------------------" <<endl
<<"Rho Matrix for Neumann BC" << endl<<
"-----------------------------------------------------------------------------"
<<endl;
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      rhoc[i][j]=0.;
      if(j == 0) rhoc[i][j] = -2.*del;
      cout<<left<<setw(8)<<rhoc[i][j];
    }
    cout<<endl;
  }
// cosine transform of rows (?)
  for (Int i=0; i<ny;i++){
      cosft1(rhoc[i]);
  }

// swap indices so we can do cosine transform on columns
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      rhotc[i][j] = rhoc[j][i];
    }
  }
  rhoc = rhotc;

// Cosine transform on columns
  for (Int i=0; i<ny;i++){
      cosft1(rhoc[i]);
  }

// swap rows back to original configuration
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      rhotc[i][j] = rhoc[j][i];
    }
  }
  rhoc = rhotc;

// solve for uh
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      uhc[i][j] = rhoc[i][j]/2./(cos(pi*i/nx)+cos(pi*j/ny)-2);
    }
  }

// sine transform uh rows
  for (Int i=0; i<ny;i++){
      cosft1(uhc[i]);
  }

// swap rows and make above the inverse sine transform
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      uc[i][j] = uhc[j][i]*2/nx;
    }
  }
  uh = u;

// sine transform uh columns
  for (Int i=0; i<ny;i++){
      cosft1(uhc[i]);
  }

// Swap rows and columns again. Make above inverse
  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      uc[i][j] = uhc[j][i]*2/nx;
    }
  }

// output the answer
cout<< endl<<
"-----------------------------------------------------------------------------" << endl
<<"Cosine Transform" << endl <<
"-----------------------------------------------------------------------------"
<<endl;

  for (Int i=0; i<nx; i++){
    for (Int j=0; j<nx; j++){
      if (abs(uc[i][j])<1e-15) uc[i][j]=0;
      cout<< left<< setw(8)<<uc[i][j];
      myfile    << std::setprecision(5)<< uc[i][j]<< " ";
    }
      cout<<endl;
      myfile <<endl;
  }

  myfile.close();


  return 0;
};
