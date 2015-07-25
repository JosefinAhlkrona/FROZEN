echo "Compiling ElmerGrid"
echo "###################"
export CC="gcc"
export CXX="g++"
export FC="gfortran"
export F77="gfortran"
export CFLAGS="-O2"
export CXXFLAGS="-O2"
export FCFLAGS="-O2"
export F90FLAGS="-O2"
export F77FLAGS="-O2"
export FFLAGS="-O2"
cd elmergrid
 make clean; ./configure --with-64bits=yes --prefix=$ELMER_HOME  && make clean; make -j4 && sudo make install
cd ..
