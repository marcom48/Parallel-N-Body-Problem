/******************************************************************************

COMP90025 Parallel and Multicore Computing Project 3 2020

Marco Marasco

This file contains helper functions for the project.

******************************************************************************/

#include "util.h"


/*
* Function used for timing program.
*/
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

/*
* Function generates a passed number of bodies within the simulation space.
* Masses (m) are initialised (in kg) with:
*       0 < m < OUT_OF_BOUNDS / 1000
* Positions (x, y, z) are initialised (in m) with:
*       -OUT_OF_BOUNDS <= x, y, z < OUT_OF_BOUNDS.
* Velocities (x, y, z) are initialised (in m/s) with 
*       -OUT_OF_BOUNDS / 100 <= x, y, z < OUT_OF_BOUNDS / 100.
*/
void generate_bodies(double total_bodies, double* body_masses, body_data *bodies){
    for(int i = 0; i < total_bodies; i++){

            // Mass of particle
            body_masses[i] = (OUT_OF_BOUNDS / 1000)*(double)rand() / (RAND_MAX + 1.0);
            
            // Compensate for zero mass.
            while (body_masses[i] == 0.0){
                body_masses[i] = (OUT_OF_BOUNDS / 1000)*(double)rand() / (RAND_MAX + 1.0);
            }
            
            // Position in x, y, z axes
            bodies[i][0] = 2 * OUT_OF_BOUNDS*(double)rand() / (RAND_MAX + 1.0) - OUT_OF_BOUNDS;
            bodies[i][1] = 2 * OUT_OF_BOUNDS*(double)rand() / (RAND_MAX + 1.0) - OUT_OF_BOUNDS;
            bodies[i][2] = 2 * OUT_OF_BOUNDS*(double)rand() / (RAND_MAX + 1.0) - OUT_OF_BOUNDS;
            
            // Velocity in x,y,z axes
            bodies[i][3] = 2 * (OUT_OF_BOUNDS / 100)*(double)rand() / (RAND_MAX + 1.0) - (OUT_OF_BOUNDS / 100);
            bodies[i][4] = 2 * (OUT_OF_BOUNDS / 100)*(double)rand() / (RAND_MAX + 1.0) - (OUT_OF_BOUNDS / 100);
            bodies[i][5] = 2 * (OUT_OF_BOUNDS / 100)*(double)rand() / (RAND_MAX + 1.0) - (OUT_OF_BOUNDS / 100);
        
        }
}

/*
* Function outputs the data of the bodies to stdout.
*/
void print_bodies(body_data *bodies, double *body_masses, double total_bodies){
    printf("\tMass\t\tx_pos\t\ty_pos\t\tz_pos\t\tx_vel\t\ty_vel\t\tz_vel\n");
    for (int i = 0; i < (int)total_bodies; i++) {
        printf("Body %d: %6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\t%6.6lf\n", i, body_masses[i],
                                        bodies[i][0], bodies[i][1], bodies[i][2],
                                        bodies[i][3], bodies[i][4], bodies[i][5]);
    }
    printf("\n\n");
}