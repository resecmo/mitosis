#include <iostream>
#include <armadillo>

using namespace arma;

//diffusion constants and matrix
const double diff_A = 1;
const double diff_Aact = 1;
const double diff_AA = 0.6;
const mat33 D = {{diff_A,0,0}, {0,diff_Aact,0}, {0,0,diff_AA}};

//ppase concentration and reaction constants
const double ppase = 1;
const double kcat = 1; //self-activation
const double kaf = 1; //aa* complex association
const double kar = 1; //aa* complex dissociation
const double kacat = 1; //aa* -> 2a* phosphorylation
const double kpf = 1; //a*ppase complex association
const double kpr = 1; //a*ppase complex dissociation
const double kpcat = 1; //appase -> a + ppase dephosphorylation

const double len = 1;
const int L = 10; //number of nodes
const double h = len / (L - 1);
vec3 conc[L];

const double tau = 0.01;

vec3 f(vec3 u) { // ppase is not accounted for
    return vec({- kcat * u[0] - kaf * u[0] * u[1],
		kcat * u[0] - kaf * u[0] * u[1] + (2 * kacat + kar) * u[2],
		kaf * u[0] * u[1] - (kacat + kar) * u[2]});
}

mat33 jac(vec3 u) {
    return mat({{-kcat - kaf * u[1], -kaf * u[0], 0},
		{kcat + kaf * u[1],  -kaf * u[0], 2 * kacat + kar},
		{kaf * u[1],         kaf * u[0],  - kacat - kar}});
}

mat33 A(vec3 u) {
    return eye(3, 3) / 12 - D * tau / (h*h) - tau / 12 * jac(u); 
}

mat33 B(vec3 u) {
    return -(eye(3, 3) * 5 / 6 + 2 * D * tau / (h*h) - tau * 5 / 6 * jac(u));
}

mat33 C(vec3 u) {
    return eye(3, 3) / 12 - D * tau / (h*h) - tau / 12 * jac(u); 
}

vec3 phi(vec3 u_prev, vec3 u_cur, vec3 u_next) {
    vec3 u_part = - u_prev / 12 - u_cur * 5 /6 - u_next / 12;
    vec3 f_part = tau * (f(u_prev) / 12 + f(u_cur) * 5 / 6 + f(u_next) / 12);
    vec3 jac_part = - tau * (jac(u_prev) * u_prev / 12 + jac(u_cur) * u_cur * 5 / 6 + jac(u_next) * u_next / 12);
    return u_part + f_part + jac_part;
} 

void initialize() {
    for (auto &vector: conc) {
        vector = vec({1,1,1});
    }
}

void step(){  
    vec3 new_conc[L];
    mat33 p[L];
    vec3 q[L];

    p[0] = eye(3, 3); // why?
    q[0] = zeros(3, 1); //why?
    for (int i = 1; i < L; i++) {
        p[i] = - inv(A(conc[i-1]) * p[i-1] - B(conc[i])) * C(conc[i+1]);
        q[i] = inv(A(conc[i-1]) * p[i-1] - B(conc[i])) * (phi(conc[i-1], conc[i], conc[i+1]) - A(conc[i-1]) * q[i-1]); // think about choosing args for phi
    }

    new_conc[L-1] = q[L-1]; // change
    for (int i = L - 2; i > -1; i--){
        new_conc[i] = p[i] * new_conc[i+1] + q[i];
    }

    std::copy(new_conc, new_conc+(L-1), conc);
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
