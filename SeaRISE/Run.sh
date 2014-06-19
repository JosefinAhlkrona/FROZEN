#!/bin/bash -l

#take so that simulation time, time step size, output interval, restart position and extruded mesh levels matches between Spinup.sif and Future.sif

#number of cores
antal=1

#partition the mesh - one piece for each core
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
echo Constructing mesh
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
ElmerGrid 2 2 mesh2D_smooth10 -partition $antal 1 1 1

# Read in geometry data for netcdf-files
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
echo Initializing Geometry Data
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
echo Initialization.sif > ELMERSOLVER_STARTINFO 
mpirun -np $antal ElmerSolver_mpi

# This is supposed to be a 125 000 ybp of only SIA, 
# initially 100 years of surface relaxation, the rest fixed surface 
# evolving temperature, forcing from GRIP using PDD, but it takes to long so for now it's just a couple of timesteps since I'm too lazy to change the files...
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
echo Paleoclimatic SpinUp
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
echo SpinUp.sif > ELMERSOLVER_STARTINFO 
mpirun -np $antal ElmerSolver_mpi

#This relaxes the surfaces a bit, its supposed to be 50 years but right now we dont want to wait for that to run but since I'm too lazy I keep the file and just to a few timesteps
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
echo Relaxing the Surface with FFS
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
echo RelaxationFFS.sif > ELMERSOLVER_STARTINFO 
mpirun -np $antal ElmerSolver_mpi

#This is the actual run, I just do 3 timesteps now so that its short
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
echo Running into the future
echo ---------------------------------------------------------------------
echo ---------------------------------------------------------------------
# 100 years into future with FFS
# climate forcing from AR4
# no temperature evolution
echo Future.sif > ELMERSOLVER_STARTINFO 
mpirun -np $antal ElmerSolver_mpi



