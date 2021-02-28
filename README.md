# Parallel N-Body Problem

> This project introduces a parallelised version of a sequential algorithm to solve the cardinal N-Body problem. The N-Body problem, published by Sir Isaac Newton in his Pricipia, considers the interaction and movements of bodies within a given space, particularly, looking at how they each influence each other’s behaviour. Solutions for N &le; 2 can be solved analytically, but for N > 2, the computational complexity becomes infeasible, andsimulation methods are required. The N-Body problem has applications in a plethora of research domains, including astrophysics, structural biology, and even machine learning. This project utilises a HPC platform to investigate the different methods to parallelise the algorithm, further using the OpenMPI and OpenMP APIs for distributed/shared memory parallelisations. The final algorithm produced by this paper achieved a 20&times; speed up using 12 nodes with 2 cores each, compared to the sequential algorithm, indicating a near linear speed up with respect to the number of cores.


## Execution
### Sequential Program
```bash
mpicxx -o sequential Sequential.cpp util.cpp -fopenmp -O3
./sequential <n_bodies> <time_steps> <secs_between_steps>
```

### Parallel Program
```bash
mpicxx -o parallel Parallel.cpp util.cpp -fopenmp -O3
mpirun -np <processes> ./parallel <n_bodies> <time_steps> <secs_between_steps>
```
