#include "/home/osetr/Desktop/mitosis/cylinder/kode/parameters.h"
#include<armadillo>
#include<iostream>

// time-independent
arma::vec::fixed<6> u_ij[I+2][J+2]; //u_ij[i][j] stands for u_i-0.5,j-0.5
//arma::vec::fixed<6> ukiter_ij[I][J]; //k-th Newton iteration
arma::sp_mat dmin_mx=arma::sp_mat();
arma::sp_mat dmax_mx=arma::sp_mat();
arma::sp_mat dav_mx=arma::sp_mat();
arma::sp_mat dzer_mx=arma::sp_mat();
arma::sp_mat C_i[I];
arma::sp_mat A_i[I];
arma::sp_mat F_i[2];
arma::sp_mat G_i[2];
arma::sp_mat BJL_ij[I][J];

// time-dependent
arma::mat::fixed<6,6> BJF_ij[I][J];
arma::mat::fixed<6,6> B_ij[I][J];
arma::mat::fixed<6,6> Binv_ij[I][J];

//initialize time-independent and initial condition
//initial condition TODO
void initializeDiffusion(){
	for (int i=0; i<3; ++i){ //should zeros be initialized???
		dav_mx(i,i)=dav;
		dmax_mx(i,i)=dmax;
		dmin_mx(i,i)=dmin;
	}
}

//1D->2D
arma::sp_mat& diffr_coeff(const int& i, const int& j){
	if (i>1 && i<=I) {return dmin_mx;}
	if (i==1) {return dav_mx;}
	return dzer_mx;
} // diffr_coeff(i,j) stands for d_i,j+0.5
arma::sp_mat& diffz_coeff(const int& i, const int& j){
	if (j>0 && j<=J){
		if (i>1) {return dmin_mx;}
		if (i==1) {return dmax_mx;}
	}
	return dzer_mx;
} // diffz_coeff(i,j) stands for d_i+0.5,j

void initializeCAFGBJL();
//1D->2D
arma::sp_mat& C_ij(const int& i, const int& j);
arma::sp_mat& A_ij(const int& i, const int& j);
arma::sp_mat& F_ij(const int& i, const int& j);
arma::sp_mat& G_ij(const int& i, const int& j);

//TODO
void newtonStep(const arma::sp_mat dt_phi[I][J]){ //arrays are passed by reference by default, I think
	int jinit=0;
	for (int i=1; i<=I; ++i){
		for (int j=jinit+1; j<=J; j+=2){ //can't proceed unless stop criterion is clear
		}
		jinit = (++jinit) % 2;
	}
	jinit=1;
	for (int i=1; i<=I; ++i){
		for (int j=jinit+1; j<=J; j+=2){
		}
		jinit = (++jinit) % 2;
	}
}
//TODO
void impEuler();
//TODO
//don't forget B_ij and Binv_ij stepwise initialization
void timeStep();

int main(){
	// larmadillo usage
	/*
	arma::vec v({1,2,3});
	arma::mat m({{0,1,0},{1,0,0},{0,0,1}});
	arma::sp_mat sm(3,3);
	sm(0,0)=0.333;
	sm(1,1)=0.5;
	sm(2,2)=1;
	std::cout << m*v << std::endl;
	std::cout << sm*v << std::endl;
	*/
	return 0;
}
