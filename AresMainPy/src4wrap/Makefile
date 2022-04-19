#!###############################################
#!#########Wrap AresMain
##########! AresMainPy
#!###############################################

#=======================================================================
#                   define the compiler names
#=======================================================================
CC       = icc
# F90      = gfortran
F90      = ifort
# F90      =  /opt/intel/composer_xe_2015.3.187/bin/intel64/ifort
PYTHON   = python
#=======================================================================
#                     additional flags
#=======================================================================
ifeq ($(F90),gfortran)
	FPP      = gfortran -E
	FPP_F90FLAGS = -x f95-cpp-input -fPIC -fpp -DMPI
	F90FLAGS = -fPIC -fpp -DMPI
    FCOMP    = gfortran
    LIBS     =
endif

ifeq ($(F90),ifort)
	FPP      = gfortran -E # gfortran f90wrap temp files only. not compilation
	FPP_F90FLAGS = -x f95-cpp-input -fPIC -cpp -DMPI
	F90FLAGS = -fpscomp logicals -fPIC -fpp -DMPI # use 1 and 0 for True and False
    FCOMP    = intelem # for f2py
	MPIFC = mpiifort
    LIBS =
endif

LIB_PATH_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/ares_lib/lib
LIB_PATH=-L${LIB_PATH_ROOT}
LIB_MKL_ROOT=${MKLROOT}/lib/intel64
LIB_MKL_PATH=-L${LIB_MKL_ROOT}
FFTW_DIR_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/ares_lib/lib/fftw
FFTW_DIR=-L${FFTW_DIR_ROOT}/lib

FFLAGS = -O3 -fPIC -fpp  -I${LIB_PATH_ROOT}/include  -module .
MPI_FFLAGS = -fPIC -O3 -traceback -heap-arrays 1000 -align array32byte -align rec32byte -fpp -DMPI -I$(FFTW_DIR_ROOT)/include -I${LIB_PATH_ROOT}/include -module .

LIBS=-lfftw3 -larpack_ifort -lxcf90 -lxc
MKLLIBS=-lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
MPI_LIBS=-larpack_ifort -lxcf90 -lxc
MPI_FFT=-lfftw3_mpi -lfftw3 
MPI_MKLLIBS=-lmkl_blacs_intelmpi_lp64 -lmkl_lapack95_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

LAPACKLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/lapack-3.10.0 -llapack -lrefblas  -ltmglib
BLASLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/BLAS-3.10.0 -lblas

CFLAGS = -fPIC #     ==> universal for ifort, gfortran, pgi

ifeq (${UNAME}, Darwin)
  LIBTOOL = libtool -static -o
else
  LIBTOOL = ar src
endif

# Python module name
PYTHON_MODN = AresMainPy
# mapping between Fortran and C types
KIND_MAP = kindMap.json

LIBSRC_SOURCES =Constants m_time_evaluate Parameters \
				Smpi_math_module array_io MathSplines Math \
				Fourier Struct_module Psp_module Read_module \
				Succeed_module Scala Grid_module Ewald \
				Finite_module XC_functional sgfft poisson_isf \
	 			Begin_module Smearing_module Nonlocalpot_module \
	 			IonLocalPotentialAssignment Isolate_module\
     			Potential_module Lapack_module Matvec_module \
				Arpack_module Chebyshev_fliter Energy_module \
	 			Mixer_module Scf_module Out_module ForceStress_module \
     			Write_module dyn_module brent cg_vasp cg_module \
     			Relax_module End_module Bands_module AresMainAPI

# file names
LIBSRC_FILES = $(addsuffix .f90,${LIBSRC_SOURCES})

# object files
LIBSRC_OBJECTS = $(addsuffix .o,${LIBSRC_SOURCES})

# only used when cleaning up
LIBSRC_FPP_FILES = $(addsuffix .fpp,${LIBSRC_SOURCES})

# names (without suffix), f90 sources
LIBSRC_WRAP_SOURCES = $(LIBSRC_SOURCES)

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90,${LIBSRC_WRAP_SOURCES})

# object files
LIBSRC_WRAP_OBJECTS = $(addsuffix .o,${LIBSRC_WRAP_SOURCES})

# fpp files
LIBSRC_WRAP_FPP_FILES = $(addsuffix .fpp,${LIBSRC_WRAP_SOURCES})

.SUFFIXES: .f90 .fpp

.PHONY: all clean

all: _${PYTHON_MODN}.so _${PYTHON_MODN}_pkg.so

clean:
	-rm -f ${LIBSRC_OBJECTS} ${LIBSRC_FPP_FILES} libsrc.a _${PYTHON_MODN}*.so \
	_${PYTHON_MODN}_pkg.so *.mod *.fpp f90wrap*.f90 f90wrap*.o *.o ${PYTHON_MODN}.py
	-rm -rf ${PYTHON_MODN}_pkg
	-rm -rf src.*/ .f2py_f2cmap .libs/ __pycache__/

.f90.o:
	#${F90} ${F90FLAGS} -c $< -o $@
	$(MPIFC) $(MPI_FFLAGS) -c $< -o $@

.c.o:
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.mod

.f90.fpp:
	${FPP} ${FPP_F90FLAGS} $<  -o $@

libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ $?

_${PYTHON_MODN}.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m ${PYTHON_MODN} ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v
	f2py-f90wrap --debug --fcompiler=$(FCOMP) --f90exec=$(MPIFC) --f77exec=$(MPIFC) \
	--f90flags='${MPI_FFLAGS}' --build-dir . -c -m _${PYTHON_MODN} f90wrap*.f90 \
	-L. -lsrc $(LAPACKLIB) $(BLASLIB) ${FFTW_DIR} ${MPI_FFT} ${LIB_PATH} $(MPI_LIBS) ${LIB_MKL_PATH} ${MPI_MKLLIBS}  


_${PYTHON_MODN}_pkg.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m ${PYTHON_MODN}_pkg ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v -P
	f2py-f90wrap --debug --fcompiler=$(FCOMP) --f90exec=$(MPIFC) --f77exec=$(MPIFC) \
	--f90flags='${MPI_FFLAGS}' --build-dir . -c -m _${PYTHON_MODN}_pkg f90wrap*.f90 \
	-L. -lsrc $(LAPACKLIB) $(BLASLIB) ${FFTW_DIR} ${MPI_FFT} ${LIB_PATH} $(MPI_LIBS) ${LIB_MKL_PATH} ${MPI_MKLLIBS} 