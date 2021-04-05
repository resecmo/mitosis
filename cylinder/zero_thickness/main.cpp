#include <iostream>
#include <armadillo>

using namespace arma;

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

const double len = 1;
const int L = 10; //number of nodes
const double h = len / (L - 1);
vec4 conc[L];

const double tau = 0.01;

vec4 f(vec4 u) { 
    return vec({- kcat*u[0] - kaf*u[0]*u[1] + kar*u[2] + kpcat*u[3],
                kcat*u[0] - kaf*u[0]*u[1] + kar*u[2] + 2*kacat*u[2] - kpf*u[1]*ppase + kpr*u[3], 
                kaf*u[0]*u[1] - kar*u[2] - kacat*u[2], 
                kpf*ppase*u[1] - kpr*u[3] - kpcat*u[3]});
}

mat44 jac(vec4 u) {
    return mat({{  -u[1]*kaf - kcat,             -u[0]*kaf,           kar,        kpcat},
	        {  -u[1]*kaf + kcat, -u[0]*kaf - ppase*kpf, 2*kacat + kar,          kpr},
	        {          u[1]*kaf,              u[0]*kaf,  -kacat - kar,            0},
	        {                 0,             ppase*kpf,             0, -kpcat - kpr}});
}

mat44 A(vec4 u) {
    return eye(4, 4) / 12 - D * tau / (h*h) - tau / 12 * jac(u); 
}

mat44 B(vec4 u) {
    return -(eye(4, 4) * 5 / 6 + 2 * D * tau / (h*h) - tau * 5 / 6 * jac(u));
}

mat44 C(vec4 u) {
    return eye(4, 4) / 12 - D * tau / (h*h) - tau / 12 * jac(u); 
}

vec4 phi(vec4 u_prev, vec4 u_cur, vec4 u_next) {
    vec4 u_part = - u_prev / 12 - u_cur * 5 /6 - u_next / 12;
    vec4 f_part = tau * (f(u_prev) / 12 + f(u_cur) * 5 / 6 + f(u_next) / 12);
    vec4 jac_part = - tau * (jac(u_prev) * u_prev / 12 + jac(u_cur) * u_cur * 5 / 6 + jac(u_next) * u_next / 12);
    return u_part + f_part + jac_part;
} 

void initialize() {
    for (auto &vector: conc) {
        vector = vec({1,1,1,1});
    }
}

void step(){  
    vec4 new_conc[L];
    mat44 p[L];
    vec4 q[L];

    p[0] = eye(4, 4); // why?
    q[0] = zeros(4, 1); //why?
    for (int i = 1; i < L; i++) {
        p[i] = - inv(A(conc[i-1]) * p[i-1] - B(conc[i])) * C(conc[i+1]);
        q[i] = inv(A(conc[i-1]) * p[i-1] - B(conc[i])) * (phi(conc[i-1], conc[i], conc[i+1]) - A(conc[i-1]) * q[i-1]); // think about choosing args for phi
    }

    new_conc[L-1] = q[L-1]; // change
    for (int i = L - 2; i > -1; i--){
        new_conc[i] = p[i] * new_conc[i+1] + q[i];
    }

    std::copy(new_conc, new_conc+L, conc);
    //return new_conc;
}

int main() {
    //std::cout << inv(mat({{3,5,1}, {0,7,0}, {9,6,1}})) << std::endl;

    initialize();
    step();    
   
    for (auto &v : conc) {
        std::cout << v << std::endl;
    }
    return 0;
}
