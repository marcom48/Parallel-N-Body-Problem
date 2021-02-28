/******************************************************************************

COMP90025 Parallel and Multicore Computing Project 3 2020

Marco Marasco

This file contains a sequential approach to solving the N-body problem
for N particles. The simulation environment is physically bounded within
a cube.

******************************************************************************/
#include "util.h"

int main(int argc, char *argv[])
{

    // Seed random.
    srand(1);

    if (argc != 4)
    {
        std::cout << "Please run file as mpirun -np <procs> ./prog <n_bodies> <simulation_step> <seconds_between_steps>." << std::endl;
        MPI_Finalize();
        return 0;
    }

    // Number of bodies for simulation.
    double total_bodies;

    // Number of logical "steps" the simulation will perform.
    double simulation_step;

    // The difference in time between each logical simulation step.
    double time_diff;

    // Holds start time of execution.
    uint64_t start;

    // Read in values.
    total_bodies = atof(argv[1]);
    simulation_step = atof(argv[2]);
    time_diff = atof(argv[3]);

    // Begin timer.
    start = GetTimeStamp();

    // Heap allocating data to allow for scale to large test cases.
    // I recognise this is slower than the stack, but the main focus for
    // analysis for this project is concerned with parallelisation techniques/APIs.

    // Maintain a separate structure to store masses.
    // Body masses are constant, thus, no need to continually broadcast
    // the data throughout the simulation.
    double *body_masses = new double[(int)total_bodies];

    // Data structure to store the bodies for the simulation.
    body_data *bodies = new body_data[(int)total_bodies];

    // Data structure to store the forces acting on a body.
    force_data *forces = new force_data[(int)total_bodies];

    // Commit a contiguous data type of 6 doubles.
    // Stores same execution related data as the body_data type.

    // Generate bodies.
    generate_bodies(total_bodies, body_masses, bodies);

    if (DEBUG)
    {
        printf("\t\t\t\t\t\tInitial body data\n");
        print_bodies(bodies, body_masses, total_bodies);
    }

    // Intialise variables required for simulation.
    int i, j;
    double x_delta, y_delta, z_delta, euclidean_distance, gravity;

    // Simulate body movements
    for (int count = 0; count < simulation_step; count++)
    {

        // Initialise force acting on each body to zero.
        // Could collapse
        for (i = 0; i < total_bodies; i++)
        {
            for (j = 0; j < 3; j++)
            {
                forces[i][j] = 0.0;
            }
        }

        // Compute forces in dimension acting on each body.
        // Could collapse this
        for (i = 0; i < total_bodies; i++)
        {
            for (j = 0; j < (int)total_bodies; j++)
            {

                // Current body working on.

                // Assert not comparing with self.
                if (j != i)
                {

                    // Displacement in x axis.
                    x_delta = bodies[i][0] - bodies[j][0];

                    // Displacement in y axis.
                    y_delta = bodies[i][1] - bodies[j][1];

                    // Displacement in z axis.
                    z_delta = bodies[i][2] - bodies[j][2];

                    // Calculate euclidean distance between bodies.
                    euclidean_distance = sqrt((x_delta * x_delta) + (y_delta * y_delta) + (z_delta * z_delta));

                    // Force of gravity / euclidean_distance between bodies.
                    gravity = GRAV_CONST * body_masses[i] * body_masses[j] / pow(euclidean_distance, 3);

                    // Compute force between bodies.
                    x_delta *= gravity;
                    y_delta *= gravity;
                    z_delta *= gravity;

                    // Add forces acting in each axis.
                    forces[i][0] += x_delta;
                    forces[i][1] += y_delta;
                    forces[i][2] += z_delta;
                }
            }
        }

        // Update forces acting on each body.
        // #pragma omp parallel for schedule(static) collapse(2)
        for (i = 0; i < total_bodies; i++)
        {

            // Iterate over axis
            for (int j = 0; j < 3; j++)
            {

                gravity = time_diff / body_masses[i];

                // Update velocity and position.
                bodies[i][j + 3] += gravity * forces[i][j];
                bodies[i][j] += time_diff * bodies[i][j + 3];

                // Body has gone out of bounds, reverse velocity in relevant axis.
                // Note:
                //      1. This is assuming elastic collisions.
                //      2. This does not impact the intended analysis for this project,
                //         it is simply a little extension to keep the system more interesting.
                if (bodies[i][j] < -OUT_OF_BOUNDS || bodies[i][j] > OUT_OF_BOUNDS)
                {
                    bodies[i][j + 3] *= -1;
                }
            }
        }
    }

    // Print bodies.
    if (DEBUG)
    {
        printf("\t\t\t\t\tBodies after %d simulation steps\n", (int)simulation_step);
        print_bodies(bodies, body_masses, total_bodies);
    }

    printf("Time: %lu us\n", (uint64_t)(GetTimeStamp() - start));

    // Clean up.
    delete body_masses;
    delete bodies;
    delete forces;

    return 0;
}