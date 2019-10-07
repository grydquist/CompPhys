# include "nr3.h"
# include "gauss_wgts.h"

template <class T>
Doub qgleg(T &func, const Doub a, const Doub b, Int N)
{
    // Weights & points
	VecDoub x(N),w(N);

    gauleg (a,b,x,w);
	Doub xm=0.5*(b+a);
	Doub xr=0.5*(b-a);
	Doub s=0;
	for (Int j=0;j<N;j++) {
		Doub dx=xr*x[j];
		s += w[j]*func(xm+dx);
	}
	return s *= xr;
}
