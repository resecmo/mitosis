#include "parameters.h"
#include<armadillo>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<string>
#include<vector>
#include<array>

// time-independent
arma::vec::fixed<6> u_ij[I+2][J+2]; //u_ij[i][j] stands for u_i-0.5,j-0.5
//arma::vec::fixed<6> ukiter_ij[I][J]; //k-th Newton iteration
arma::sp_mat dbulk_mx=arma::sp_mat(6,6);
arma::sp_mat dmt_mx=arma::sp_mat(6,6);
arma::sp_mat dav_mx=arma::sp_mat(6,6);
arma::sp_mat dzer_mx=arma::sp_mat(6,6);
arma::sp_mat C_i[I];//these 4 matrices will have function(i,j) interfaces
arma::sp_mat A_i[I];
arma::sp_mat F_i[3];
arma::sp_mat G_i[3];
//jacobianless
arma::sp_mat BJL_ij[I][J];
arma::sp_mat E_mx(6,6);

// time-dependent
// jacobianful
arma::mat::fixed<6,6> BJF_ij[I][J];
arma::mat::fixed<6,6> B_ij[I][J];
arma::mat::fixed<6,6> Binv_ij[I][J];

//initialize time-independent and initial condition
auto astest(const int& i, const int& j){//toRemove
	return u_ij[i+1][j+1][1];
}
void set_init_cond_from_file(const std::string& path){ //JxI.csv
	std::ifstream istream(path);
	std::string line;
	int i=0;
	int j=0;
	int k=0;
	if (istream.is_open()){
		while (getline(istream,line)){
			if (j>=J){throw "J value exceeded";}
			std::stringstream ss(line);
			std::vector<std::vector<double>> us;
			//std::string sus;
			std::vector<double> u;
			std::string su;
			std::string _su;
			while (getline(ss, su, ';')){
				if (i>=I){throw "I value exceeded";}
				std::stringstream _ss(su);
				while (getline(_ss, _su, ',')){
					if (k>=6){throw "vector of dimension gt 6";}
					u.push_back(stod(_su));
					++k;
				}
				arma::vec _u(u);
				u_ij[i+1][j+1] = _u;
				//astest = _u(1); //toRemove
				//us.push_back(u);
				u.clear();
				k=0;
				++i;
			}
			/*
			for (arma::vec _u : us){
				u_ij[i+1][j+1] = _u;
				//astest = _u(1); //toRemove
			}
			*/
			i=0;
			++j;
		}
	}
	}
void initializeDiffusion(){
	for (int i=0; i<3; ++i){ //should zeros be initialized???
		dav_mx(i,i)=dav;
		dmt_mx(i,i)=dmt;
		dbulk_mx(i,i)=dbulk;
	}
}

//1D->2D
bool r_in_mt(const int& i, const int& j){
	return i < rmti && i > 0;
}
bool r_in_bulk(const int& i, const int& j){
	return i > rmti && i < I; //was <=
}
bool r_in_av(const int& i, const int& j){//is incorrect and not called
	return i == rmti;
}
bool r_in_vanish(const int& i, const int& j){
	return i==0 || i==I; // was I+1
}
arma::sp_mat& diffr_coeff(const int& i, const int& j){
	if (r_in_bulk(i,j)) {return dbulk_mx;}
	if (r_in_mt(i,j)) {return dmt_mx;}
	if (r_in_vanish(i,j)) {return dzer_mx;}
	return dav_mx;
} // diffr_coeff(i,j) stands for d_i,j+0.5
bool z_in_mt(const int& i, const int& j){
	return i < rmti && j > 0 && j < J; //was <=
}
bool z_in_bulk(const int& i, const int& j){
	return i >= rmti && j > 0 && j < J; //was <=
}
bool z_in_av(const int& i, const int& j);
bool z_in_vanish(const int& i, const int& j){
	return j == 0 || j == J; // was J+1
}
arma::sp_mat& diffz_coeff(const int& i, const int& j){//??? 
	if (z_in_bulk(i,j)) {return dbulk_mx;}
	if (z_in_mt(i,j)) {return dmt_mx;}
	if (z_in_vanish(i,j)) {return dzer_mx;}
	return dav_mx;
	/*
	if (j > 0 && j <= J){
		if (i >= rmti) {return dbulk_mx;}
		return dmt_mx;
	}
	return dzer_mx;
	*/
} // diffz_coeff(i,j) stands for d_i+0.5,j

void initializeCAFGEBJL(){
	for (int i=0; i<I; ++i){
		C_i[i] = (ht*(i + 1.) / (hr*hr*(i + 0.5))) * diffr_coeff(i+1,0);
		A_i[i] = (ht*i / (hr*hr*(i+0.5))) * diffr_coeff(i,0);
	}
	F_i[0] = dzer_mx;
	F_i[1] =(ht / (hz*hz)) * dbulk_mx;
	F_i[2] =(ht / (hz*hz)) * dmt_mx;
	G_i[0] = dzer_mx;
	G_i[1] =(ht / (hz*hz)) * dbulk_mx;
	G_i[2] =(ht / (hz*hz)) * dmt_mx;
	for (int i=0; i<6; ++i){
		E_mx(i,i)=1;
	}
	for (int i=0; i<I; ++i){
		for (int j=0; j<J; ++j){
			BJL_ij[i][j] = E_mx + (ht / (hr*hr*(i+0.5))) * ((i + 1.)*diffr_coeff(i+1,j) + i*diffr_coeff(i,j))
				       + (ht / (hz*hz)) * (diffz_coeff(i,j+1) + diffz_coeff(i,j));
		}
	}

}
//1D->2D
arma::sp_mat& C_ij(const int& i, const int& j){
	return C_i[i];
}
arma::sp_mat& A_ij(const int& i, const int& j){
	return A_i[i];
}
arma::sp_mat F_ij(const int& i, const int& j){
	/* question is what return type should be (to & or not to &)?
	 * another question is who optimizes better -- me or compiler?
	if (j < J-1){
		if (i > rmti) {
			return F_i[1];
		}
		return F_i[2];
	}
	return F_i[0];
	*/
	return (j < J-1)*((i >= rmti) ? F_i[1] : F_i[2]) + dzer_mx; //was i>
}
arma::sp_mat G_ij(const int& i, const int& j){
	/* same questions here
	if (j > 0){
		if (i > rmti){
			return G_i[1];
		}
		return G_i[2];
	}
	return G_i[0];
	*/
	return (j > 0)*((i >= rmti) ? G_i[1] : G_i[2]) + dzer_mx; //was i>
}

//TODO
bool satisfied=false;
void update(const int& i, const int& j,
		arma::vec::fixed<6>& uij_prev,
		double& shift, double& worstshift,
		const arma::vec::fixed<6> dt_phi[I][J]){ //arrays are passed by reference by default, I think
	uij_prev=u_ij[i+1][j+1];
	//u_ij[i+1][j+1] = Binv_ij[i][j]
		   //* (dt_phi[i][j]
		      //- altomg*BJL_ij[i][j]*u_ij[i+1][j+1] //sign?
		      //+ C_ij(i,j)*u_ij[i+2][j+1]
		      //+ A_ij(i,j)*u_ij[i][j+1]
		      //+ F_ij(i,j)*u_ij[i+1][j+2]
		      //+ G_ij(i,j)*u_ij[i+1][j]);
	u_ij[i+1][j+1] = (1 - omega)*u_ij[i+1][j+1] 
		+ omega*Binv_ij[i][j]
		   * (dt_phi[i][j]
		      + C_ij(i,j)*u_ij[i+2][j+1]
		      + A_ij(i,j)*u_ij[i][j+1]
		      + F_ij(i,j)*u_ij[i+1][j+2]
		      + G_ij(i,j)*u_ij[i+1][j]);
	shift = norm(u_ij[i+1][j+1] - uij_prev, "inf")
		/ std::max(1., norm(uij_prev, "inf"));
	worstshift = std::max(worstshift, shift);
}
arma::vec::fixed<6> dt_phi[I][J];
void newtonStep(){ //arrays are passed by reference by default, I think
	int jinit=0;
	double worstshift=0;
	double shift=0;
	arma::vec::fixed<6> uij_prev;
	/*
	auto update = [&](const int& i, const int& j){
		uij_prev=u_ij[i+1][j+1];
		u_ij[i+1][j+1]=Binv_ij[i][j]
			   * (dt_phi[i][j] 
			      + C_ij(i,j)*u_ij[i+2][j+1]
			      + A_ij(i,j)*u_ij[i][j+1]
			      + F_ij(i,j)*u_ij[i+1][j+2]
			      + G_ij(i,j)*u_ij[i+1][j]);
		shift = norm(u_ij[i+1][j+1] - uij_prev, "inf")
			/ std::max(1., norm(uij_prev, "inf"));
		worstshift = std::max(worstshift, shift);
	};
	*/
	for (int i=0; i<I; ++i){
		for (int j=jinit; j<J; j+=2){
			update(i, j, uij_prev, shift, worstshift, dt_phi);
		}
		jinit = (jinit + 1) % 2;
	}
	jinit=1;
	for (int i=0; i<I; ++i){
		for (int j=jinit; j<J; j+=2){
			update(i, j, uij_prev, shift, worstshift, dt_phi);
		}
		jinit = (jinit + 1) % 2;
	}
	satisfied = worstshift < nepsilon;
}
//TODO
int max_newton_steps=0;
int newton_steps=0;
arma::vec::fixed<6> u11; //toRemove
double _u11;
void impEuler(){
	for (int i=0; i<I; ++i){
		for (int j=0; j<J; ++j){
			const arma::vec::fixed<6> uij = u_ij[i+1][j+1];
			BJF_ij[i][j] = - ht*ourJac(uij,i,j); //mind the sign
			//B_ij[i][j] = invomg*BJL_ij[i][j] + BJF_ij[i][j];//omega or invomg???
			B_ij[i][j] = BJL_ij[i][j] + BJF_ij[i][j];
			Binv_ij[i][j] = B_ij[i][j].i();
			dt_phi[i][j] = ht*rhs(uij,i,j) + BJF_ij[i][j]*uij + uij; //SIGN!!!
			//u_ij[i+1][j+1] = uij; // I forgot what I meant
		}
	}
	newton_steps=0;
	do{
		newtonStep();
		++newton_steps;
		u11 = u_ij[1][1]; //toRemove
		_u11 = u11[1];
	} while (!satisfied);
	//satisfied=false;
	max_newton_steps = (newton_steps > max_newton_steps)?newton_steps:max_newton_steps;
	//std::cout << "newton_steps: " << newton_steps << std::endl;
}
//TODO
//don't forget B_ij and Binv_ij stepwise initialization
void stamp(std::ostream& csv_stream){ //was ofstream)))
	for (int j=1; j <= J; ++j){
		for (int i=1; i<I; ++i){
			csv_stream << u_ij[i][j][1] << ",";
		}
		csv_stream << u_ij[I][j][1] << std::endl;
	}
}
void timeStep(std::ostream& csv_stream){ // i wonder if ofstream should better be reopened at each call
	stamp(csv_stream);
	impEuler();
}
void saveinfo(const std::string& save_path, int argc, char* argv[]){
	std::ofstream infostream(save_path + ".info");
	/*
	std::vector<std::string> args;
	std::for_each(std::begin(argv),std::end(argv),
			[&](const std::string& word){args.push_back(word);});
	*/
	std::string word;
	if (infostream.is_open()){
		for (int i=0; i<argc-1; ++i){
			word = argv[i];
			infostream << argv[argc-1] << " ";
		}
		infostream << word << std::endl;
	} else {std::cout << "Can't save " << save_path << ".info" << std::endl;}
}

class NullBuffer : public std::streambuf{
	public:
		int overflow(int c) {return c;}
}; //stackoverflow)))
int main(int argc, char* argv[]){
	std::string help_msg = std::string("`-ic init_cond.csv` for initial condition\n")
				+std::string("`-t time` for time of simulation\n")
				+std::string("`-o save_path` for main output file\n")
				+std::string("-o option creates save_path.csv and save_path.info\n")
				+std::string("w/ u_ijs and other information respectively")
				+std::string("\n`-omg omega` for omega\n")
				+std::string("`-eps relax_crit` for nepsilon\n")
				+std::string("`-ht ht` for ht\n")
				+std::string("`-e each` for each\n");
	 
	std::string init_cond_path="cylinder_ic.scv";
	double time=10;
	int each = 1;
	std::string save_path="out";
	std::string comment="none";
	for (int i=1; i<argc; ++i){
		std::string strarg(argv[i]);
		if (strarg=="-ic"){
			init_cond_path=argv[i+1];
		}
		if (strarg=="-t"){
			time=std::stod(argv[i+1]);
		}
		if (strarg=="-o"){
			save_path=argv[i+1];
		}
		if (strarg=="-m"){ //m for message
			comment=argv[i+1];
		}
		if (strarg=="-e"){
			each=std::stod(argv[i+1]);
		}
		if (strarg=="-omg"){
			omega=std::stod(argv[i+1]);
		}
		if (strarg=="-ht"){
			ht=std::stod(argv[i+1]);
		}
		if (strarg=="-eps"){
			nepsilon=std::stod(argv[i+1]);
		}
		if (strarg=="-h"){
			std::cout << help_msg;
			return 0;
		}
	}
	unsigned long long timesteps = time/ht;
	initializeDiffusion();
	initializeCAFGEBJL();
	set_init_cond_from_file(init_cond_path);
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
	//actual execution finally starts
	std::ofstream csv_stream(save_path + ".csv");
	NullBuffer nullBuf;
	std::ostream nullStr(&nullBuf); //stackoverflow)))
	if (csv_stream.is_open()){
		unsigned long long t=0;
		for (t=0; t<timesteps; ++t){
			if (t % each == 0)
				timeStep(csv_stream);
			else
				timeStep(nullStr); //stackoverflow)))
		}
		if (t % each == 0)
			stamp(csv_stream);
		csv_stream.close();
	} else {std::cout << "Can't save " << save_path << ".csv" << std::endl;}
	saveinfo(save_path, argc, argv);
	std::cout << "max Newton steps: " << max_newton_steps << std::endl;
	std::cout << "number of timesteps: " << timesteps << std::endl;
	return 0;
}
