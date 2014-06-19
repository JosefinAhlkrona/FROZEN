#!/bin/bash -l

elmerf90 -o SIASolverJosefin2 LinearSliding.f90 SIASolverJosefin2.f90
elmerf90 -o FlowSolveSIAFS ErrorEstimationSubs.f90 FlowSolveSIAFS.f90
elmerf90 -o linkZsGreen_newData linkZsGreen_newData.f90
elmerf90 -o LinearSliding LinearSliding.f90

sudo mv  SIASolverJosefin2 /usr/local/Elmer/lib/
sudo mv  FlowSolveSIAFS /usr/local/Elmer/lib/
sudo mv  linkZsGreen_newData /usr/local/Elmer/lib/
sudo mv  LinearSliding /usr/local/Elmer/lib/

