/******************************************************************************

COMP90025 Parallel and Multicore Computing Project 3 2020

Marco Marasco

This file contains a distributed parallelised approach to solving the
N-body problem for N particles, using OpenMP and OpenMPI. The simulation
environment is physically bounded within a sphere

The approach is to distribute jobs to each worker, who will compute the
results, and then return them to the root node to compute the final values
required for output.

File variance: Reduces message size by splitting position and velocity data.

******************************************************************************/
#include "utilVel.h"


// Initialise communicator constant.
const MPI_Comm comm = MPI_COMM_WORLD;



int main(int argc, char *argv[]) {

    // Seed random.
    srand(1);

    // Initialise MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Number of bodies for simulation.
    double total_bodies;

    // Number of logical "steps" the simulation will perform.
    double simulation_step;

    // The difference in time between each logical simulation step.
    double time_diff;

    // Send/recv buffer for first broadcast.
    double start_data[3];
    
    // Holds start time of execution.
    uint64_t start;

    if (rank == root){

        // Read in values.
        total_bodies = atof(argv[1]);
        simulation_step = atof(argv[2]);
        time_diff = atof(argv[3]);


        // Begin timer.
        start = GetTimeStamp();

        // Broadcast the number of bodies, simulation steps, 
        // and time difference to all workers.
        start_data[0] = total_bodies;
        start_data[1] = simulation_step;
        start_data[2] = time_diff;
        MPI_Bcast(start_data, 3, MPI_DOUBLE, 0, comm);
    
    } else{
        
        // Receive the number of bodies, simulation steps, 
        // and time difference from the master.
        MPI_Bcast(start_data, 3, MPI_DOUBLE, 0, comm);

        total_bodies = start_data[0];
        simulation_step = start_data[1];
        time_diff = start_data[2];
    }

    if ((int)total_bodies % size != 0){
            if (rank == root)
                printf("Please ensure that the number of bodies is divisible by the number of MPI workers.\n");
            MPI_Finalize();
            return 0;
    }
    
   
    // Heap allocating data to allow for scale to large test cases.
    // I recognise this is slower than the stack, but the main focus for
    // analysis for this project is concerned with parallelisation techniques/APIs.


    // Maintain a separate structure to store masses.
    // Body masses are constant, thus, no need to continually broadcast
    // the data throughout the simulation.    
    double *body_masses = new double[(int)total_bodies];
    
    // Data structure to store the body positions for the simulation.
    body_data *positions = new body_data[(int)total_bodies];

    // Data structure to store the body positions for the simulation.
    body_data *velocities = new body_data[(int)total_bodies];

    // Number of bodies processor will handle.
    int processor_bodies = (int)total_bodies / size;

    // Create a pointer to the first body the procedure will compute.
    // Improves readability by removing complex index arithmetic.
    body_data *proc_pos = positions + rank * processor_bodies;
    body_data *proc_vel = velocities + rank * processor_bodies;

    // Data structure to store the forces acting on a body.
    force_data *forces =  new force_data[processor_bodies]; 
    
  
    // Commit a contiguous data type of 6 doubles.
    // Stores same execution related data as the body_data type.
    MPI_Datatype MPI_body_data;
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_body_data);
    MPI_Type_commit(&MPI_body_data);

    if (rank == root) {
        
        // Generate bodies.
        generate_bodies(total_bodies, body_masses, positions, velocities);

        if (DEBUG){
            printf("\t\t\t\t\t\tInitial body data\n");
            print_bodies(positions, velocities, body_masses, total_bodies);

        }

    }

    
    // Brodcast the masses of the bodies.
    MPI_Bcast(body_masses, total_bodies, MPI_DOUBLE, 0, comm);
    
    // Broadcast the initial data for the bodies.
    MPI_Bcast(positions, total_bodies, MPI_body_data, 0, comm);
    MPI_Bcast(velocities, total_bodies, MPI_body_data, 0, comm);
    


    // Intialise variables required for simulation.
    int i, j, curr_body;
    double x_delta, y_delta, z_delta, euclidean_distance, gravity;

    // Simulate body movements
    for (int count = 0; count < simulation_step; count++) {
        

        // Initialise force acting on each body to zero.
        // Could collapse
        for(int i = 0; i < processor_bodies; i++){
            forces[i][0] = 0.0;
            forces[i][1] = 0.0;
            forces[i][2] = 0.0;
        }
        

        // Compute forces in dimension acting on each body.
        // Could collapse this
        for (i = 0; i < processor_bodies; i++) {
                for (j = 0; j < (int)total_bodies; j++) {
                    
                    // Current body working on.
                    curr_body = rank * processor_bodies + i;

                    // Assert not comparing with self.
                    if (j != curr_body) {
                    
                        // Displacement in x axis.
                        x_delta = proc_pos[i][0] - positions[j][0];

                        // Displacement in y axis.
                        y_delta = proc_pos[i][1] - positions[j][1];

                        // Displacement in z axis.
                        z_delta = proc_pos[i][2] - positions[j][2];

                        // Calculate euclidean distance between bodies. 
                        euclidean_distance = sqrt((x_delta * x_delta) + (y_delta * y_delta) + (z_delta * z_delta));
                        
                        // Force of gravity / euclidean_distance between bodies.
                        gravity = GRAV_CONST * body_masses[curr_body] * body_masses[j] / pow(euclidean_distance, 3);

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
        for (i = 0; i < processor_bodies; i++) {

            // Iterate over axis
            for(int j = 0; j < 3; j++){

                curr_body = rank * processor_bodies + i;
                
                gravity = time_diff / body_masses[curr_body];

                // Update position and velocity.
                proc_pos[i][j] += time_diff * proc_vel[i][j];
                proc_vel[i][j] += gravity * forces[i][j];

                // Body has gone out of bounds, reverse velocity in relevant axis.
                // Note:
                //      1. This is assuming elastic collisions.
                //      2. This does not impact the intended analysis for this project,
                //         it is simply a little extension to keep the system more interesting.
                if (proc_pos[i][j] < -OUT_OF_BOUNDS || proc_pos[i][j] > OUT_OF_BOUNDS){
                    proc_vel[i][j] *= -1;
                }
                
            }

            

        }

        // Gather positions and velocities.
        MPI_Allgather(MPI_IN_PLACE, processor_bodies, MPI_body_data, positions, processor_bodies, MPI_body_data, comm);
        MPI_Allgather(MPI_IN_PLACE, processor_bodies, MPI_body_data, velocities, processor_bodies, MPI_body_data, comm);

        // Print bodies.
        // if (DEBUG && rank == root) {
        //     printf("\t\t\t\t\tBodies at simulation step %d/%d\n", count + 1, (int)simulation_step);
        //     print_bodies(bodies, body_masses, total_bodies);
               
        // }


    }

    // Print bodies.
    if (rank == root) {

        if (DEBUG){
            printf("\t\t\t\t\tBodies after %d simulation steps\n", (int)simulation_step);
            print_bodies(positions, velocities, body_masses, total_bodies);
        }

        printf ("Time: %lu us\n", (uint64_t) (GetTimeStamp() - start));   
    }

    // Clean up.
    delete body_masses;
    delete positions;
    delete velocities;
    delete forces;
    MPI_Type_free(&MPI_body_data);
    MPI_Finalize();
    
    return 0;
}