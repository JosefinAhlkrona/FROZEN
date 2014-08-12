#!/bin/bash -f
#

#mkdir trunk

# this works for Debian/Ubuntu/LinuxMint
export MPIDIR="/usr/lib/openmpi"


export CC="mpicc"
export CXX="mpicxx"
export FC="mpif90"
export F77="mpif90"
export F90="mpif90"

#export CC=gcc
#export CXX=g++
#export FC=gfortran
#export F77=gfortran
#export F90="gfortran"
#export CPP=g++

# replace with preferred installation directory

#ELMER_REV="rev$(svn info trunk|grep Revision|awk '{print $2}')"
export ELMER_ROOT="/usr/local"
#export ELMER_INSTALL="${ELMER_ROOT}/${ELMER_REV}"
export ELMER_INSTALL="/usr/local/ElmerRev6822"

# Trilinos installation location

trilinos_install="/usr/include/trilinos-11.2.5/"


echo "ELMER installation directory:" 
echo $ELMER_INSTALL


# flags are for GNU compiler
# use -03, if you want to push performance (and your luck)
export OPTFLAGS=""
export CFLAGS="$OPTFLAGS"
export CXXFLAGS="$OPTFLAGS"
export FCFLAGS="$OPTFLAGS"
export F90FLAGS="$OPTFLAGS"
export F77FLAGS="$OPTFLAGS"
export FFLAGS="$OPTFLAGS"


export MUMPS="/usr"

# CHOLMOD stuff (that ent broken with the latest Ubuntu update - but not needed)
#export CHOLMOD_INC="/usr/include/suitesparse/"
#export CFLAGS="$CFLAGS -DHAVE_CHOLMOD -I$CHOLMOD_INC"
#export FCPPFLAGS="$FCPPFLAGS -DHAVE_CHOLMOD"
#export LDFLAGS="$LDFLAGS -L/usr/lib -lcholmod -lcamd -lccolamd -lcolamd -lamd -lmetis"

# TRILINOS stuff (the painfull part) proofed to work with libtrilinos
# 10.4.0.dfsg-1ubuntu2-dev
export CXXFLAGS="$CXXFLAGS -DHAVE_TRILINOS  -I${trilinos_install}lib -I${trilinos_install}include"

export FCPPFLAGS="$FCPPFLAGS -DHAVE_TRILINOS"
#export LDFLAGS="$LDFLAGS -L/usr/lib -ltrilinos_belostpetra -ltrilinos_belosepetra -ltrilinos_belos -ltrilinos_ml -ltrilinos_ifpack -ltrilinos_amesos -ltrilinos_galeri -ltrilinos_isorropia -ltrilinos_epetraext -ltrilinos_tpetrainout -ltrilinos_tpetra -ltrilinos_triutils -ltrilinos_zoltan -ltrilinos_epetra -ltrilinos_kokkoslinalg -ltrilinos_kokkosnodeapi -ltrilinos_kokkos -ltrilinos_teuchos"
export LDFLAGS="$LDFLAGS -L/usr/lib -L${trilinos_install}lib -lbelostpetra -lbelosepetra -lbelos -lml -lifpack -lamesos -lgaleri -lisorropia -lepetraext -ltpetrainout -ltpetra -ltriutils -lzoltan -lepetra -lkokkoslinalg -lkokkosnodeapi -lkokkos -lteuchosremainder -lteuchoscomm -lteuchoscore -lteuchosnumerics -lteuchosparameterlist"




modules="matc umfpack mathlibs meshgen2d eio hutiter fem"



#cd trunk
##### configure, build and install #########
 for m in $modules; do
   echo "compiling module $m"
   echo "#############################"
   cd $m
   pwd
      make clean; ./configure   --with-64bits=yes --with-mpi=yes --with-hypre="-lHYPRE" --with-mumps="-I$MUMPS/include -L$MUMPS/lib -ldmumps_ptscotch -lmumps_common_ptscotch" --with-lapack="-L/usr/lib -llapack" --with-blas="-L/usr/lib -lblas" --prefix=$ELMER_INSTALL  
#      make clean; ./configure    --with-hypre="-lHYPRE" --with-mumps="-I$MUMPS/include -L$MUMPS/lib -ldmumps_ptscotch -lmumps_common_ptscotch" --with-lapack="-L/usr/lib -llapack" --with-blas="-L/usr/lib -lblas" --prefix=$ELMER_INSTALL 
	make -j4  && echo "Installing into $ELMER_INSTALL";  make install
   cd ..
done

