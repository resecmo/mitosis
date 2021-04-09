#include "parameters.h"
#include <iostream>
#include <armadillo>

using namespace arma;


const double len = 1;
const int L = 10; //number of nodes
const double h = len / (L - 1);
vec4 conc[L];

const double tau = 0.01;


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
    vec4 u_part = u_prev / 12 + u_cur * 5 /6 + u_next / 12; 
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

    p[0] = eye(4, 4);
    q[0] = zeros(4, 1);
    for (int i = 1; i < L-1; i++) {
        p[i] = - inv(A(conc[i-1]) * p[i-1] - B(conc[i])) * C(conc[i+1]);
        q[i] = inv(A(conc[i-1]) * p[i-1] - B(conc[i])) * (phi(conc[i-1], conc[i], conc[i+1]) - A(conc[i-1]) * q[i-1]);
    }

    // A_last * u[L-2] + B_last * u[L-1] == phi_last
    mat44 A_last = eye(4, 4) / 6 - 2 * tau * D / (h*h) - tau * jac(conc[L-2]) / 6;
    mat44 B_last = eye(4, 4) * 5 / 6 + 2 * tau * D / (h*h) - tau * jac(conc[L-1]) * 5 / 6;
    vec4 phi_last = conc[L-2] / 6 + conc[L-1] * 5 / 6 + tau * (f(conc[L-2]) / 6 + f(conc[L-1]) * 5 / 6 - jac(conc[L-2]) * conc[L-2] / 6 - jac(conc[L-1]) * conc[L-1] * 5 / 6);
    new_conc[L-1] = inv(A_last * p[L-2] + B_last) * (phi_last - A_last * q[L-2]);
    for (int i = L - 2; i > -1; i--){
        new_conc[i] = p[i] * new_conc[i+1] + q[i];
    }

    std::copy(new_conc, new_conc+L, conc);
}

int main() {

    initialize();
    step();    
    
    for (auto &v : conc) {
        std::cout << v << std::endl;
    }

    return 0;
}
