#include <armadillo>
#include <cmath>

static double nepsilon=1e-05; //epsilon for newton iterations
const int I=100;
const int J=100;
const double R=1;
const double Z=1;
const double RMT=0;
const double hr=R/I;
const double hz=Z/J;
static double ht=1e-02;
const int rmti=I*(RMT/R); //might require +-1 correction
static double omega = 1;
const double invomg = 1/omega;
//const double altomg = 1 - invomg;
const double altomg = omega - 1;

const double kaf=0.042;
const double kar=1;
const double kac=0.0014;
const double km=0.001;//???
const double kpf=0.013;
const double kpr=0.025;
const double kpc=0.0013;
const double knf=25.3;
const double knr=842;
const double kpos=0.001;//???
const double kneg=0.001;//???
const double knegst=0.001;//???

const double alpha=1;
//const double n0=1;
double fgus(const int& i, const int& j, const double& sg=10, const double& mur=0, const double& muz=50){
	//return std::exp(-0.5 * (std::pow(mur - i, 2) + std::pow(muz - j, 2)) / std::pow(sg, 2));
	return 1;
}
double n0(const int& i, const int& j){
	return fgus(i,j);
}

const double dbulk=0.1;
const double dmt=0.01;
const double dav=2 * dbulk * dmt / (dbulk+dmt);

//following are to optimize jacobian calculation
const double j12_ = -1 - kpc*km;
//const double knfn0 = knf*n0;
double knfn0(const int& i, const int& j){
	return knf * n0(i,j);
}
const double knfal0 = knf*alpha;
//const double j22_ = kpr*km - kpf - knfn0;
double j22_(const int& i, const int& j){
	return kpr*km - kpf - knfn0(i,j);
}
const double j23_ = kar + 2*kac;
const double j242 = knfal0 - kpos;
const double j26_ = kneg + knegst;
const double j33_ = - (kar + kac);
const double j56_ = kneg + knegst;
const double j66_ = - (2*kneg + kpos);
//jacobian
//arma::mat::fixed<6,6> ourJac(const arma::vec::fixed<6>& u){
	//arma::mat res = {
		//{0,0,0,0,0,0},
		//{0,0,0,0,0,0},
		//{0,0,0,0,0,0},
		//{0,0,0,0,0,0},
		//{0,0,0,0,0,0},
		//{0,0,0,0,0,0}
	//};
	//return res;
//}

//double _knfan456(const arma::vec::fixed<6>& u, const int& i, const int& j){
	//return knfn0(i,j) - knfal0*(u[3] + u[4] + u[5]);
//}
arma::mat::fixed<6,6> ourJac(const arma::vec::fixed<6>& u, const int& i, const int& j){
	const double knfan456 = knfn0(i,j) - knfal0*(u[3] + u[4] + u[5]);
	//const double knfan456 = _knfan456(u,i,j);
	const double knfau1 = knfal0*u[0];
	const double knfau2 = knfal0*u[1];
	const double kafu1 = kaf*u[0];
	const double kafu2 = kaf*u[1];
	arma::mat res = {
		{-1 - knfan456 - kaf*u[1] - kpos*u[4], j12_, kar, knr - knfau1, knfau1, kneg + knfau1},
		{1 - kafu2, j22_(i,j) - kafu1 - knfan456 - kpos*u[2], j23_, j242*u[1], knr + knfau2, j26_ + knfau2},
		{kafu2,		    kafu1,		 j33_,	       0,	 0,	   0},
		{knfan456,    -kpos*u[3],    0,    -(knr + kpos*u[1] + knfau1),	 -knfau1,    kneg-knfau1},
		{-kpos*u[4],    knfan456,    0,	 -knfau2,   -(knr + kpos*u[0] + knfau2),    j56_ - knfau2},
		{kpos*u[4],   kpos*u[3],    0,	   kpos*u[1],	 kpos*u[0],	 j66_}
	};
	return res;
}

//right hand side (chemical reactions)
//arma::vec::fixed<6> rhs(const arma::vec::fixed<6>& u){
	//arma::vec::fixed<6> res = {0,0,0,0,0,0};
	//return res;
//}

arma::vec::fixed<6> rhs(const arma::vec::fixed<6>& u, const int& i, const int& j){
	const double kafu12=kaf*u[0]*u[1];
	const double knfan456 = knfn0(i,j) - knfal0*(u[3] + u[4] + u[5]);
	arma::vec::fixed<6> res = {
		(-1 - kpos*u[4] - knfan456)*u[0] - kafu12 + kpc*km*u[1] + kar*u[2] + knr*u[3] + kneg*u[5],
		u[0] - kafu12 + (kpr*km - kpf - kpos*u[3] - knfan456)*u[1] + 2*kac*u[2] + knr*u[4] + (kneg + knegst)*u[5],
		kafu12 - (kar + kac)*u[2],
		knfan456*u[0] - (knr + kpos*u[1])*u[3] + kneg*u[5],
		knfan456*u[1] - (knr + kpos*u[0])*u[4] + (kneg + knegst)*u[5],
		-(2*kneg + kpos)*u[5] + kpos*(u[0]*u[4] + u[1]*u[3])
	};
	return res;
}

