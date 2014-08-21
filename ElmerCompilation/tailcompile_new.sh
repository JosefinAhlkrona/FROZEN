#! /bin/bash

##### re-link #####
echo "re-linking Elmer directory"
echo "##########################"
export ELMER_ROOT="/usr/local"
export ELMER_REV="ElmerRev6822"

pushd $ELMER_ROOT
pwd
echo "linking $ELMER_REV"
sudo chmod -R a+rX $ELMER_REV
sudo rm -f Elmer
sudo ln -s $ELMER_REV Elmer
ls -ltrh
export ELMER_HOME="$ELMER_ROOT/Elmer"
export PATH="$PATH:$ELMER_HOME/bin"
echo "target directory: $ELMER_HOME"
popd

##### compile Elmer/Ice ######
echo "Compiling Elmer/Ice"
echo "###################"
cd elmerice
export ELMER_ROOT="/usr/local"
export ELMER_HOME=${ELMER_ROOT}"/Elmer"
make purge
make compile
sudo -E make install
cd netcdf2
elmerf90 GridDataReader.f90 -o GridDataReader.so
sudo mv GridDataReader.so $ELMER_HOME/share/elmersolver/lib/
sudo chmod a+rx,g+rwx $ELMER_HOME/share/elmersolver/lib/GridDataReader.so
cd ../..
pwd
# ParStokes
echo "Compiling ParStokes"
echo "###################"
cd fem/src/modules/
pwd
cp ParStokes.src ParStokes.f90
elmerf90 ParStokes.f90 -o ParStokes.so
sudo mv ParStokes.so $ELMER_HOME/share/elmersolver/lib/
sudo chmod a+rx,g+rwx $ELMER_HOME/share/elmersolver/lib/ParStokes.so
cd ../../..
##### compile ElmerGrid ###### !USE compile_grid.sh instead
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

export CC="mpicc"
export CXX="mpicxx"
export FC="mpif90"
export F77="mpif90"
export F90="mpif90"

cd elmergrid
 ./configure --with-64bits=yes --prefix=$ELMER_HOME  && make clean; make -j4 && sudo make install
cd ..
##### compile ElmerPost ####### And this one does work instead of compile_post.sh
echo "Compiling ElmerPost"
echo "###################"
export CFLAGS="-O2"
export CXXFLAGS="-O2"
export FCFLAGS="-O2"
export F90FLAGS="-O2"
export F77FLAGS="-O2"
export FFLAGS="-O2"
#cd trunk
cd post

# ./configure --with-64bits=yes --prefix=$ELMER_HOME  && make clean; make -j4 && sudo make install

cd src/plugins
chmod u+x ./configure
export CFLAGS="-I/usr/share/ffmpeg -I/usr/include/tcl8.5 -I/usr/include/libavcodec/ -I/usr/include/libswscale"
export LIBS="-lavcodec -lavutil -lswscale"
export ELMER_POST_HOME=$ELMER_HOME/share/elmerpost

 ./configure --prefix=$ELMER_POST_HOME && make && sudo make install

cd ../../..
##### compile ElmerGUI
echo "Compiling ElmerGUI"
echo "###################"

echo ""
pwd
echo ""

#cp ElmerGUI.pri trunk/ElmerGUI
cd ElmerGUI
make clean
qmake
make -j4
sudo make install
cd ../..
echo "ALL DONE"
