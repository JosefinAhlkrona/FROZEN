#ElmerGrid 14 2 cirkel.msh -autoclean #convert gmsh grid into a elmer type grid

#ElmerSolver partitiontest.sif #get a division into FS and SIA nodes

ElmerSolver findweights.sif    #convert this info into nodal weights

ElmerGrid 2 2 cirkel -metis 4 0  #partionion the mesh according to the weights

echo lookatmesh.sif > ELMERSOLVER_STARTINFO #run something to get to see the partioning in paraview
mpirun -np 4 ElmerSolver_mpi
