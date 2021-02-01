#include "/home/osetr/Desktop/mitosis/kode/parameters.h"
#include<armadillo>
#include<iostream>

// time-independent
arma::vec::fixed<6> u_ij[I+2][J+2]; //u_ij[i][j] stands for u_i-0.5,j-0.5
arma::sp_mat dr_ij[I+1][J+2]; //dr_ij stands for d_i,j-0.5, flow along r
arma::sp_mat dz_ij[I+2][J+1]; //dz_ij stands for d_i-0.5,j, flow along z
arma::mat::fixed<6,6> C_ij[I][J];
arma::mat::fixed<6,6> A_ij[I][J];
arma::mat::fixed<6,6> F_ij[I][J];
arma::mat::fixed<6,6> G_ij[I][J];
arma::mat::fixed<6,6> BJL_ij[I][J];

// time-dependent
arma::mat::fixed<6,6> BJF_ij[I][J];
arma::mat::fixed<6,6> B_ij[I][J];

//TODO
void initializeInitial(); //initialize time-independent and inital condition
//TODO
void newtonStep();
//TODO
void impEuler();
//TODO
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
