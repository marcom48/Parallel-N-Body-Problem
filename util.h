/******************************************************************************

COMP90025 Parallel and Multicore Computing Project 3 2020

Marco Marasco

This file contains helper items for the project.

******************************************************************************/
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

/*
body_data is a helper type definition, it stores (respectively):
x position, y position, z position, x velocity, y velocity, z velocity.

It has been implemented to improve readability, of the code, removing
more complex index arithmetic to access body data.
*/
typedef double body_data[6];


/*
force_data is a helper type definition, it stores (respectively):
x axis force, y axis force, z axis force for a given body.

It has been implemented to improve readability, of the code, removing
more complex index arithmetic to access force data.
*/
typedef double force_data[3];

/*
Gravitational constant for computing force. This is relatively arbitrary with respect
to the intended outcomes of the analysis this project. The unit is m^(3) * kg^(-1) * s^(-2).
*/
const double GRAV_CONST = 60.0;
// const double GRAV_CONST = 6.67408e-11; // Actual gravitational constant

// Max distance from an axis a body can me (metres).
#define OUT_OF_BOUNDS 1000.0

// Debugging flag.
#define DEBUG false

// Flag for root worker in MPI communicator.
const int root = 0;


// Timer function.
uint64_t GetTimeStamp();

// Data generation function.
void generate_bodies(double total_bodies, double* body_masses, body_data *postions);

// Debugging function.
void print_bodies(body_data *bodies, double *body_masses, double total_bodies);