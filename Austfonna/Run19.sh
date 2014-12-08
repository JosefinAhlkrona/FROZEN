#!/bin/bash -l

#SBATCH -A p2013035
#SBATCH -p core -n 16
#SBATCH -t 08:00:00
#SBATCH -J Austfonna_test

cores=1



module load elmer

rm mesh2D/init.result.*

#elmerf90 GridDataReader.f90 -o GridDataReader.so  -lnetcdf  -lnetcdff -I $ELMER_HOME/share/netcdf/include
#elmerf90 -o linkZsGreen_newData linkZsGreen_newData.f90

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/comp/intel/Compiler/composer_xe_2013.3.163/mkl/lib/intel64

ElmerGrid 2 2  mesh2D -partition ${cores} 1 1 1 -autoclean

echo initialization19.sif > ELMERSOLVER_STARTINFO 


mpirun -np ${cores} ElmerSolver_mpi

wait

echo austfonna19.sif > ELMERSOLVER_STARTINFO 


mpirun -np ${cores} ElmerSolver_mpi

