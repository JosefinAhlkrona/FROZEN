#!/bin/bash -l

#debug=" -g -fbounds-check -mcmodel=large -fbacktrace"
debug="-O2"

elmerf90 -o Functionals ${debug}   Functionals.f90
elmerf90 -o ErrorEstimationSubs ${debug}  ErrorEstimationSubs.f90
elmerf90 -o NavierStokes2 ${debug}   NavierStokes2.f90
elmerf90 -o LinearSliding ${debug} LinearSliding.f90
elmerf90 -o SIASolverJosefin2 ${debug}   LinearSliding.f90 SIASolverJosefin2.f90
##elmerf90 -o FlowSolveSIAFS Functionals.f90 ErrorEstimationSubs.f90 FlowSolveSIAFS.f90
elmerf90 -o FlowSolveSIAFS ${debug}   Functionals.f90 ErrorEstimationSubs.f90 NavierStokes2.f90 FlowSolveSIAFS.f90
##elmerf90 -o FlowSolveSIAFS2 Functionals.f90 ErrorEstimationSubs.f90 FlowSolveSIAFS2.f90
##mpif90 ${debug} -O0 -m64 -fPIC -fPIC -I. -Ibinio -I../binio -I/usr/include -m64 -fPIC -I/usr/local/ElmerRev6697/share/elmersolver/include -shared -o FlowSolveSIAFS Functionals.f90 ErrorEstimationSubs.f90 FlowSolveSIAFS.f90 
elmerf90 -o ComputeNormal ${debug}  ComputeNormal.f90
elmerf90 -o linkZsGreen_newData ${debug} linkZsGreen_newData.f90
elmerf90 -o  NeighbourFinder ${debug}  NeighbourFinder.f90 

#debug=" -g -fbounds-check -mcmodel=large -fbacktrace"
elmerf90 -o Grid2DInterpolator2 ${debug}  Grid2DInterpolator.f90

elmerfolder=ElmerRev7138

sudo mv -f SIASolverJosefin2 /usr/local/${elmerfolder}/lib/
sudo mv  -f FlowSolveSIAFS /usr/local/${elmerfolder}/lib/
#sudo mv  -f FlowSolveSIAFS2 /usr/local/${elmerfolder}/lib/
sudo mv -f  ComputeNormal /usr/local/${elmerfolder}/lib/
sudo mv -f  linkZsGreen_newData /usr/local/${elmerfolder}/lib/
sudo mv -f  LinearSliding /usr/local/${elmerfolder}/lib/

sudo mv -f  NeighbourFinder /usr/local/${elmerfolder}/lib/

sudo mv -f  Grid2DInterpolator2 /usr/local/${elmerfolder}/lib/
