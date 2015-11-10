#!/bin/bash -l

#SBATCH -A m.2015-1-311

#SBATCH -J frozen

#SBATCH -t 1:00:00

##SBATCH --ntasks-per-node=32

#SBATCH -N 1

#SBATCH -n 16

#SBATCH -e error_file.e 

#SBATCH -o output_file.o 

. /opt/modules/default/init/bash

module add cmake/2.8.12.2

module add cce/8.3.4

module add cray-libsci fftw

module add cray-trilinos/11.10.1.0

module add cray-hdf5/1.8.13

module add cray-netcdf/4.3.2

module add  cray-tpsl/1.4.2

module swap PrgEnv-cray PrgEnv-gnu

ulimit -s unlimited 

echo

ntask=`env | grep SLURM | grep SLURM_NTASKS`

echo

cores=`echo ${ntask} | awk -F"=" '{print $2}'`


echo "the number of cores are ${cores}"

export CRAY_ROOTFS=DSL

first_letter=`echo ${USER} | cut -c1`

export ELMER_HOME1="/cfs/nobackup/""${first_letter}""/""${USER}""/elmer"

export ELMER_HOME="/cfs/nobackup/s/saefa/FROZEN/elmer"



echo ${ELMER_HOME}



#export ELMER_HOME="/cfs/nobackup/""${first_letter}""${USER}""/FROZEN/elmer"
export PATH=$ELMER_HOME/bin:$ELMER_HOME/lib:$ELMER_HOME/include:$ELMER_HOME1/bin:$ELMER_HOME1/lib:$ELMER_HOME1/include:$PATH
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ELMER_HOME/lib:$ELMER_HOME/include:"/cfs/nobackup/s/saefa/FROZEN/elmerfem/fem/src":"/opt/cray/tpsl/1.4.2/GNU/49/x86_64/lib":"/usr/lib64"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ELMER_HOME1/lib:$ELMER_HOME1/include:$ELMER_HOME/lib:$ELMER_HOME/include:"/opt/cray/tpsl/1.4.2/GNU/49/x86_64/lib":"/usr/lib64":"/opt/intel/composer_xe_2013_sp1.4.211/compiler/lib/intel64"

export LIBS="-L$ELMER_HOME/lib -L$ELMER_HOME1/lib -L/usr/lib64":"/opt/intel/composer_xe_2013_sp1.4.211/compiler/lib/intel64"
#export LDFLAGS="-L$ELMER_HOME/lib -L/cfs/nobackup/s/saefa/FROZEN/elmerfem/fem/src -L/usr/lib64"

export LDFLAGS="-L$ELMER_HOME/lib -L$ELMER_HOME1/lib -L/usr/lib64"

##cores="12"

#rm mesh2D/test15fs.*

#elmerf90 GridDataReader.f90 -o GridDataReader.so  -lnetcdf  -lnetcdff -I $ELMER_HOME/share/netcdf/include
#elmerf90 -o linkZsGreen_newData linkZsGreen_newData.f90

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/comp/intel/Compiler/composer_xe_2013.3.163/mkl/lib/intel64

#strace ElmerGrid 2 2 ./mesh2D -partition ${cores} 1 1 1 

#-autoclean

echo ${PATH}

echo initialization25.sif > ELMERSOLVER_STARTINFO 


aprun -n ${cores} ElmerSolver_mpi   

wait


echo test25fsISCAL.sif > ELMERSOLVER_STARTINFO 


#mpirun -np ${cores} valgrind --leak-check=full --track-origins=yes --leak-check=full --show-leak-kinds=all  ElmerSolver_mpi


aprun -n ${cores} ElmerSolver_mpi 


#> my_output_file 2>&1



#mpirun -np ${cores} ElmerSolver_mpi


