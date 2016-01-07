#!/bin/bash -l

. /opt/modules/default/init/bash

module add cmake/2.8.12.2

module add cce/8.3.4

#module add PrgEnv-cray/5.2.40                                                                                                                

#module add PrgEnv-gnu/5.2.40                                                                                                                 


#module add cray-petsc                                                                                                                        

module add cray-libsci fftw

module add cray-trilinos/11.10.1.0

module add cray-tpsl/1.5.0                                                                                                                   

module add cray-hdf5/1.8.13

#module add cray-hdf5-parallel/1.8.13                                                                                                         

module add cray-netcdf/4.3.2




module swap PrgEnv-cray PrgEnv-gnu

#export CRAYPE_LINK_TYPE=dynamic                                                                                                              





#/opt/cray/trilinos/default/CRAY/83/x86_64/lib                                                                                                



#######################!/bin/sh -f                                                                                                            


export CRAYPE_LINK_TYPE=dynamic


export CRAYPE_LINK_TYPE=dynamic


export ELMER_HOME="/cfs/nobackup/n/ninakir/elmer"
export PATH=$ELMER_HOME/bin:$ELMER_HOME/lib:$ELMER_HOME/include:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ELMER_HOME/lib:$ELMER_HOME/include:"/cfs/nobackup/n/ninakir/elmerfem/fem/src":"/opt/cray/tpsl/1\
.5.0/GNU/49/x86_64/lib":"/opt/intel/composer_xe_2013_sp1.4.211/compiler/lib/intel64"
export LIBS="-L$ELMER_HOME/lib"
export LDFLAGS="-L$ELMER_HOME/lib -L/cfs/nobackup/n/ninakir/elmerfem/fem/src -L/usr/lib64 -L/opt/intel/composer_xe_2013_sp1.4.211/compiler/lib/intel64"
#export LDFLAGS="-R /cfs/nobackup/s/saefa/FROZEN/elmerfem/fem/src/"                                                                           


export CC="cc"
export CXX="CC"
export FC="ftn"
export F77="ftn"
export F90="ftn"


elmerf90 -o Functionals Functionals.f90
elmerf90 -o ErrorEstimationSubs ErrorEstimationSubs.f90
elmerf90 -o NavierStokes2 NavierStokes2.f90
elmerf90 -o LinearSliding LinearSliding.f90
elmerf90 -o SIASolverJosefin2 LinearSliding.f90 SIASolverJosefin2.f90
#elmerf90 -o FlowSolveSIAFS Functionals.f90 ErrorEstimationSubs.f90 FlowSolveSIAFS.f90
elmerf90 -o FlowSolveSIAFS Functionals.f90 ErrorEstimationSubs.f90 NavierStokes2.f90 FlowSolveSIAFS.f90
#elmerf90 -o FlowSolveSIAFS2 Functionals.f90 ErrorEstimationSubs.f90 FlowSolveSIAFS2.f90
#mpif90 -g -O0 -m64 -fPIC -fPIC -I. -Ibinio -I../binio -I/usr/include -m64 -fPIC -I/usr/local/ElmerRev6697/share/elmersolver/include -shared -o FlowSolveSIAFS Functionals.f90 ErrorEstimationSubs.f90 FlowSolveSIAFS.f90 
elmerf90 -o ComputeNormal ComputeNormal.f90
elmerf90 -o linkZsGreen_newData linkZsGreen_newData.f90
elmerf90 -o Grid2DInterpolator2   Grid2DInterpolator.f90

firstletter="n"
username="ninakir"

mv  SIASolverJosefin2 /cfs/nobackup/${firstletter}/${username}/elmer/lib
mv  FlowSolveSIAFS /cfs/nobackup/${firstletter}/${username}/elmer/lib
mv  FlowSolveSIAFS2 /cfs/nobackup/${firstletter}/${username}/elmer/lib
mv  ComputeNormal /cfs/nobackup/${firstletter}/${username}/elmer/lib
mv  linkZsGreen_newData /cfs/nobackup/${firstletter}/${username}/elmer/lib
mv  LinearSliding /cfs/nobackup/${firstletter}/${username}/elmer/lib
mv Grid2DInterpolator2 /cfs/nobackup/${firstletter}/${username}/elmer/lib

rm -f *.mod

