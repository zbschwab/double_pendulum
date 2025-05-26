#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "phys_math.h"

#define MAXTIME 60  // in seconds
#define T_STEP 0.01  // timestep dt in seconds
#define PI 3.1415927
#define G 9.8

// mechanics precision conventions: rad needs 6 digits, kg and meter need 3

// converts angle from degrees to radians
double deg_to_rad(double theta_deg) {
    double theta_rad;
    theta_rad = theta_deg * PI / 180;
    return theta_rad;
}

// converts angle from polar to cartesian coordinates
void polar_to_car(double theta, double length, cart_t *coords) {
    coords->x = length * sin(theta);
    coords->y = length * cos(theta);
}

// takes state and returns the derivative of each state variable
deriv_t deriv(state_t *s, constants_t *c) {
    // calculate f and alpha for coupled ODEs (see equation sheet)
    double f_1, f_2, alpha_1, alpha_2, denom, M, diff;

    M = (c->m_2 / (c->m_1 + c->m_2));
    diff = s->theta_1 - s->theta_2;

    f_1 = - (c->l_2 / c->l_1) * M * pow(s->omega_2, 2) * sin(diff) - (G / c->l_1 * sin(s->theta_1));
    f_2 = (c->l_1 / c->l_2) * pow(s->omega_1, 2) * sin(diff) - (G / c->l_2 * sin(s->theta_2));

    alpha_1 = (c->l_2 / c->l_1) * M * cos(diff);
    alpha_2 = (c->l_1 / c->l_2) * cos(diff);
    denom = 1 - (alpha_1 * alpha_2);

    deriv_t d = {0};
    d.dtheta_1 = s->omega_1;
    d.dtheta_2 = s->omega_2;
    d.d2theta_1 = (f_1 - alpha_1 * f_2) / denom;
    d.d2theta_2 = (- alpha_2 * f_1 + f_2) / denom;

    return d;
}

// integrates using RK4 method to predict the state for given timestep
void runge_kutta_4(state_t *s, constants_t *c, double t) {
    state_t temp = *s;

    // k1 = f(t, y)
    deriv_t k1;
    k1 = deriv(s, c);

    // k2 = f(t + dt/2, temp)
    // k2 = f(t + dt/2, y + dt/2 * k1)
    deriv_t k2;
    temp.theta_1 = s->theta_1 + 0.5*T_STEP*k1.dtheta_1;
    temp.theta_2 = s->theta_2 + 0.5*T_STEP*k1.dtheta_2;
    temp.omega_1 = s->omega_1 + 0.5*T_STEP*k1.d2theta_1;
    temp.omega_2 = s->omega_2 + 0.5*T_STEP*k1.d2theta_2;
    k2 = deriv(&temp, c);
    

    // k3 = f(t + dt/2, temp)
    // k3 = f(t + dt/2, y + dt/2 * k2)
    deriv_t k3;
    temp.theta_1 = s->theta_1 + 0.5*T_STEP*k2.dtheta_1;
    temp.theta_2 = s->theta_2 + 0.5*T_STEP*k2.dtheta_2;
    temp.omega_1 = s->omega_1 + 0.5*T_STEP*k2.d2theta_1;
    temp.omega_2 = s->omega_2 + 0.5*T_STEP*k2.d2theta_2;
    k3 = deriv(&temp, c);

    // k4 = f(t + dt, temp)
    // k4 = f(t + dt, y + dt/2 * k3)
    deriv_t k4;
    temp.theta_1 = s->theta_1 + 0.5*T_STEP*k3.dtheta_1;
    temp.theta_2 = s->theta_2 + 0.5*T_STEP*k3.dtheta_2;
    temp.omega_1 = s->omega_1 + 0.5*T_STEP*k3.d2theta_1;
    temp.omega_2 = s->omega_2 + 0.5*T_STEP*k3.d2theta_2;
    k4 = deriv(&temp, c);

    // y_next = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    // t_next = t + dt (handled in main loop)
    s->theta_1 += T_STEP/6 * (k1.dtheta_1 + 2*k2.dtheta_1 + 2*k3.dtheta_1 + k4.dtheta_1);
    s->theta_2 += T_STEP/6 * (k1.dtheta_2 + 2*k2.dtheta_2 + 2*k3.dtheta_2 + k4.dtheta_2);
    s->omega_1 += T_STEP/6 * (k1.d2theta_1 + 2*k2.d2theta_1 + 2*k3.d2theta_1 + k4.d2theta_1);
    s->omega_2 += T_STEP/6 * (k1.d2theta_2 + 2*k2.d2theta_2 + 2*k3.d2theta_2 + k4.d2theta_2);
}

