#! /bin/bash

# Written by Evan Gowan

# download Trilinos from the official website: http://trilinos.org/download/

# Place this script in the same directory that contains the Trilinos source code

# before running this script, make sure you have installed (which are not by default in Linux Mint/Ubuntu):
# - cmake
# - g++
# - gfortran
# - mpi-default-dev
# - libblas-dev
# - liblapack-dev

# the file names this installs are different than the old Debian package. You have to change the buildelmer.sh script for that to work (I will upload that when I have figured out how to compile Elmer)

build_nodes=6 # number of processors you want to use to compile this

source_path=$( pwd ) # note, if you run this script from a directory other than where the source code code is located, you will have to change this

echo ${source_path}

rm -r -f builddir # in case something screwed up

mkdir builddir
cd builddir


# I put the location of the trilinos binaries in a folder called /usr/include/trilinos/
# You have to use a modified version of buildelmer.sh to be able to use this, as the file naming scheme is
# a bit different than the old Debian binary

install_directory="/usr/include/trilinos-11.2.5/"

# this is a tested serial compilation

#c_compiler=$( which gcc )
#gpp_compiler=$( which g++ )
#fortran_compiler=$( which gfortran )

#cmake \
#-DCMAKE_C_COMPILER=${c_compiler} \
#-DCMAKE_CXX_COMPILER=${gpp_compiler} \
#-DCMAKE_Fortran_COMPILER=${fortran_compiler} \
#-DTrilinos_ENABLE_ALL_PACKAGES=ON \
#-DCMAKE_INSTALL_PREFIX=/usr/include/trilinos/ \
#${source_path}



# If you want to use MPI, use this invocation instead

c_compiler=$( which mpicc )
gpp_compiler=$( which mpicxx )
fortran_compiler=$( which mpif90 )

cmake -DCMAKE_C_COMPILER=${c_compiler} -DCMAKE_CXX_COMPILER=${gpp_compiler} -DCMAKE_Fortran_COMPILER=${fortran_compiler} -DTrilinos_ENABLE_ALL_PACKAGES=ON -DCMAKE_INSTALL_PREFIX=${install_directory} -DBUILD_SHARED_LIBS:BOOL=ON -DTPL_ENABLE_MPI=ON   "${source_path}"


make -j${build_nodes} install
