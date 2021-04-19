#include <armadillo>

using namespace arma;

const double len = 1;
const int L = 10; //number of nodes
const double h = len / (L - 1);
vec4 conc[L];

const double tau = 0.01;

//diffusion constants and matrix
const double diff_A = 1;
const double diff_Aact = 1;
const double diff_AA = 0.6;
const double diff_APpase = 1;
const mat44 D = diagmat(vec({diff_A, diff_Aact, diff_AA, diff_APpase})); 

//ppase concentration and reaction constants
const double ppase = 1;
const double kcat = 1; //self-activation
const double kaf = 1; //aa* complex association
const double kar = 1; //aa* complex dissociation
const double kacat = 1; //aa* -> 2a* phosphorylation
const double kpf = 1; //a*ppase complex association
const double kpr = 1; //a*ppase complex dissociation
const double kpcat = 1; //a*ppase -> a + ppase dephosphorylation

vec4 f0(vec4 u) { 
    return vec({- kcat*u[0] - kaf*u[0]*u[1] + kar*u[2] + kpcat*u[3],
                kcat*u[0] - kaf*u[0]*u[1] + kar*u[2] + 2*kacat*u[2] - kpf*u[1]*ppase + kpr*u[3], 
                kaf*u[0]*u[1] - kar*u[2] - kacat*u[2], 
                kpf*ppase*u[1] - kpr*u[3] - kpcat*u[3]});
}

mat44 jac0(vec4 u) {
    return mat({{  -u[1]*kaf - kcat,             -u[0]*kaf,           kar,        kpcat},
	        {  -u[1]*kaf + kcat, -u[0]*kaf - ppase*kpf, 2*kacat + kar,          kpr},
	        {          u[1]*kaf,              u[0]*kaf,  -kacat - kar,            0},
	        {                 0,             ppase*kpf,             0, -kpcat - kpr}});
}

void initialize0(){
    for (int i = 0; i < L; i++) {
        conc[i] = ones<vec>(4);
    }
}



//test1, no reactions, uniform initial concentration distribution
vec4 f_test1(vec4 u) { 
    return zeros<vec>(4);
}

mat44 jac_test1(vec4 u) {
    return zeros<mat>(4, 4);
}

void initialize_test1(){
    for (int i = 0; i < L; i++) {
        conc[i] = ones<vec>(4);
    }
}

//arrays of test functions
vec4(*(f_test[]))(vec4) = {f_test1};
mat44(*(jac_test[]))(vec4) = {jac_test1};
void(*(initialize_test[]))() = {initialize_test1};

