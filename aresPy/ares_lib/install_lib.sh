#!/bin/bash

#> AUTHOR: lantian xue
#> DATE: 2021.12.29
#> EMAIL: xlt@calypso.cn


#> func: build_dir
build_dir(){
  if  [ -d $1 ] ;then
    echo "rm $1 and rebuild"
    rm -r $1 
    mkdir $1
  else
    echo "Build directory $1"
    mkdir $1
  fi
}


#> func: install ARPACK LIBXC FFTW
install_ARPACK(){
  build_dir 'ARPACK'
  tar -xzvf ARPACK.tar.gz -C ARPACK
  cd ARPACK/ARPACK

  cat ./ARmake.inc | awk '{if(NR==31){$0="home='$(pwd)'"}}\
  /^FC/{$3="'$FC'"}\
  /^PFC/{$3="'$MPIFC'"}\
  {print}' > aa

  mv aa ARmake.inc
  make clean
  make all
  cd ../..
}


install_libxc(){
  build_dir 'LIBXC'
  tar -xzvf libxc-3.0.1.tar.gz -C LIBXC
  cd ./LIBXC/libxc*
  ./configure --prefix=`pwd` --with-pic CC=$CC FC=$FC
  make
  make install
  cd ../..
}


install_fftw() {
  build_dir 'FFTW'
  tar -xzvf fftw-3.3.9.tar.gz -C FFTW
  cd FFTW/fftw*
  ./configure --prefix=`pwd` --enable-mpi --with-pic CC=$CC FC=$FC MPICC=$MPICC MPIFC=$MPIFC F77=$F77
  make
  make install
  cd ../..
}


#> copy librarys to lib
copy_libs(){
  build_dir 'lib'
  cp ARPACK/ARPACK/lib*.a ./lib

  cp LIBXC/libxc*/lib/*.a ./lib
  cp -r LIBXC/libxc*/include ./lib

  build_dir 'lib/fftw'
  cp -r FFTW/fftw*/lib ./lib/fftw
  cp -r FFTW/fftw*/include ./lib/fftw
}  


#> ------------------------Main part----------------------------
#> environmet variable
FC=ifort
MPIFC=mpiifort
CC=icc
MPICC=mpiicc
F77=ifort


echo "install ARPACK"
install_ARPACK > log_arpack
echo "ARPACK installed"

echo "install LIBXC"
install_libxc > log_xc
echo "LIBXC installed"

echo "install FFTW"
install_fftw > log_fftw
echo "FFTW installed"

echo 'copy libraties'
copy_libs
echo 'libraries installed ok'
