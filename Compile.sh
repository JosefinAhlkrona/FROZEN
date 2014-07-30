#!/bin/bash -l

elmerf90 -o SIASolverJosefin2 LinearSliding.f90 SIASolverJosefin2.f90
elmerf90 -o FlowSolveSIAFS Functionals.f90 ErrorEstimationSubs.f90 FlowSolveSIAFS.f90
#mpif90 -g -O0 -m64 -fPIC -fPIC -I. -Ibinio -I../binio -I/usr/include -m64 -fPIC -I/usr/local/ElmerRev6697/share/elmersolver/include -shared -o FlowSolveSIAFS ErrorEstimationSubs.f90 FlowSolveSIAFS.f90 
elmerf90 -o ComputeNormal ComputeNormal.f90
elmerf90 -o linkZsGreen_newData linkZsGreen_newData.f90
elmerf90 -o LinearSliding LinearSliding.f90

sudo mv  SIASolverJosefin2 /usr/local/Elmer/lib/
sudo mv  FlowSolveSIAFS /usr/local/Elmer/lib/
sudo mv  ComputeNormal /usr/local/Elmer/lib/
sudo mv  linkZsGreen_newData /usr/local/Elmer/lib/
sudo mv  LinearSliding /usr/local/Elmer/lib/


