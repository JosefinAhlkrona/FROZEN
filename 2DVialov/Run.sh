#!/bin/bash -l

antal=2

ElmerGrid 2 2 square -partition $antal 1 1

echo FFS.sif > ELMERSOLVER_STARTINFO 
mpirun -np $antal ElmerSolver_mpi



