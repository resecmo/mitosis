#include <armadillo>

const int I=100;
const int J=100;
const double R=1;
const double Z=1;
const double hr=R/I;
const double hz=Z/J;
const double ht=0.1;

const double kaf=1;
const double kar=1;
const double kac=1;
const double km=1;
const double kpf=1;
const double kpr=1;
const double kpc=1;
const double knf=1;
const double knr=1;
const double kpos=1;
const double kneg=1;
const double knegst=1;

const double alpha=1;
const double n0=1;

const double dmin=1;
const double dmax=10;
const double dav=0.5*(dmin+dmax)/(dmin*dmax);

//following are to optimize jacobian calculation
const double j12_ = -1 - kpc*km;
const double knfn0 = knf*n0;
const double knfal0 = knf*alpha;
const double j22_ = kpr*km - kpf - knfn0;
const double j23_ = kar + 2*kac;
const double j242 = knfal0 - kpos;
const double j26_ = kneg + knegst;
const double j33_ = - (kar + kac);
const double j56_ = kneg + knegst;
const double j66_ = - (2*kneg + kpos);

arma::mat::fixed<6,6> ourJac(const arma::vec::fixed<6>& u){
	const double knfan456 = knfn0 - knfal0*(u[3] + u[4] + u[5]);
	const double knfau1 = knfal0*u[0];
	const double knfau2 = knfal0*u[1];
	const double kafu1 = kaf*u[0];
	const double kafu2 = kaf*u[1];
	arma::mat res = {
		{-1 - knfan456 - kaf*u[1] - kpos*u[4], j12_, kar, knr - knfau1, knfau1, kneg + knfau1},
		{1 - kafu2, j22_ - kafu1 - knfan456 - kpos*u[2], j23_, j242*u[1], knr + knfau2, j26_ + knfau2},
		{kafu2,		    kafu1,		 j33_,	       0,	 0,	   0},
		{knfan456,    -kpos*u[3],    0,    -(knr + kpos*u[1] + knfau1),	 -knfau1,    kneg-knfau1},
		{knfan456,    -kpos*u[3],    0,	 -knfau2,   -(knr + kpos*u[0] + knfau2),    j56_ - knfau1},
		{kpos*u[4],   kpos*u[3],    0,	   kpos*u[1],	 kpos*u[0],	 j66_}
			};
	return res;
}
