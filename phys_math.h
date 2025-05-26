#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

# pragma once

typedef enum ics {  // ics = initial conditions
    M1, M2,
    L1, L2,
    T1, T2
} ics_t;

// holds user-given constants
typedef struct constants {
    double m_1;  // mass at end of 1st pendulum in kg
    double m_2;  // mass at end of 2nd pendulum in kg
    double l_1;  // length of 1st pendulum in meters
    double l_2;  // length of 2nd pendulum in meters
} constants_t;

// holds pendulum angles and their 1st derivatives
typedef struct state {
    double theta_1;
    double theta_2;
    double omega_1;
    double omega_2;
} state_t;

// holds 1st and 2nd derivatives of pendulum angles
typedef struct deriv {
    double dtheta_1;
    double dtheta_2;
    double d2theta_1;
    double d2theta_2;
} deriv_t;

// position in cartesian coordinates
typedef struct cart {
    double x;
    double y;
} cart_t;


// converts angle from degrees to radians
double deg_to_rad(double theta_deg);

// converts angle from polar to cartesian coordinates
void polar_to_car(double theta, double length, cart_t *coords);

// takes state and returns the derivative of each state variable
deriv_t deriv(state_t *s, constants_t *c);

// integrates using RK4 method to predict the state for given timestep
void runge_kutta_4(state_t *s, constants_t *c, double t);
