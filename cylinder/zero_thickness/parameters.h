#include <armadillo>

using namespace arma;

double len = 20;
int L = 400; //number of nodes
double h = len / (L - 1);
vec4* conc = 0;

double tau = 0.1;

//diffusion constants and matrix
const double diff_A = 1;
const double diff_Aact = 1;
const double diff_AA = 0.6;
const double diff_APpase = 1;
mat44 D0 = diagmat(vec({diff_A, diff_Aact, diff_AA, diff_APpase})); //TODO make const?

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

const mat44 D_test1 = diagmat(ones<vec>(4));

//test2, no reactions, initial concentration p
vec4 f_test2(vec4 u) { 
    return zeros<vec>(4);
}

mat44 jac_test2(vec4 u) {
    return zeros<mat>(4, 4);
}

void initialize_test2(){
    conc[0] = ones<vec>(4);
    for (int i = 1; i < L; i++) {
        conc[i] = zeros<vec>(4) * i;
    }
}

const mat44 D_test2 = diagmat(ones<vec>(4));

//test3, Kolmogorov-Petrovskiy-Piskunov
const double kpp_const = 2.0;
vec4 f_test3(vec4 u) { 
    return vec({kpp_const * u[0] * (1- u[0]),
                0, 
                0, 
                0});
}

mat44 jac_test3(vec4 u) {
    return mat({{ kpp_const * (1 - 2 * u[0]), 0, 0, 0},
                {                          0, 0, 0, 0},
                {                          0, 0, 0, 0},
                {                          0, 0, 0, 0}});
}

void initialize_test3(){
    conc[0] = ones<vec>(4);
    for (int i = 1; i < L; i++) {
        conc[i] = zeros<vec>(4) / 2;
    }
}

const mat44 D_test3 = diagmat(vec({0.001, 0, 0, 0}));

//test4, FitzHugh-Nagumo
const double fhn_i = 0.0;
const double fhn_a = -0.7;
const double fhn_b = 0.8;
const double fhn_phi = 0.04;
vec4 f_test4(vec4 u) { 
    return vec({u[0] - u[0]*u[0]*u[0] / 3 - u[1] + fhn_i,
                fhn_phi * (u[0] - fhn_a - fhn_b * u[1]), 
                0, 
                0});
}

mat44 jac_test4(vec4 u) {
    return mat({{ 1 - u[0] * u[0],                -1, 0, 0},
                {         fhn_phi, - fhn_b * fhn_phi, 0, 0},
                {               0,                 0, 0, 0},
                {               0,                 0, 0, 0}});
}

void initialize_test4(){
    int initial_length = (int)(0.1 / h) + 1;  // length of spike in initial conditions
    for (int i = 0; i < initial_length; i++) {
        conc[i] = vec({1.5, 0.6, 0, 0});
    }
    for (int i = initial_length; i < L; i++) {
        conc[i] = vec({-1.1994, -0.6243, 0, 0});
    }
}

const mat44 D_test4 = diagmat(vec({0.001, 0, 0, 0}));

//arrays of test functions
vec4(*(f_test[]))(vec4) = {f_test1, f_test2, f_test3, f_test4};
mat44(*(jac_test[]))(vec4) = {jac_test1, jac_test2, jac_test3, jac_test4};
void(*(initialize_test[]))() = {initialize_test1, initialize_test2, initialize_test3, initialize_test4};
const mat44* D_test[] = {&D_test1, &D_test2, &D_test3, &D_test4};

