.SUFFIXES: .f90 .o
#---------------------------------------------------------------------------------
#> (0)easy config
#---------------------------------------------------------------------------------
#> 'intel' or 'gcc'
#PLATFORM_ENV='gcc'
PLATFORM_ENV='intel'
#---------------------------------------------------------------------------------
#> (1)perpare for different platforms
#---------------------------------------------------------------------------------
#> gfortran
ifeq ($(PLATFORM_ENV),'gcc')
FC = gfortran
LINKER = $(FC)
MPIFC = mpif90
MPILINKER = $(MPIFC)

FFLAGS = -O3 -cpp  -ffree-line-length-none  -I${LIB_PATH_ROOT}/include -J${WORKDIR}
FFLAGS_debug = -g -cpp  -ffree-line-length-none  -I${LIB_PATH_ROOT}/include -J${WORKDIR}
#FFLAGS = ${FFLAGS_debug}
OMP_FFLAGS = -g -cpp  -ffree-line-length-none  -Wtabs -fbounds-check -fcheck=all -DMPI -I${FFTW_DIR_ROOT}/include -I${LIB_PATH_ROOT}/include -J${OMP_WORKDIR}
#-fsanitize=address 
#OMP_FFLAGS = -g -O3 -cpp  -ffree-line-length-none  -fcheck=all -DMPI -I${LIB_PATH_ROOT}/include -I${FFTW_DIR_ROOT}/include -J${OMP_WORKDIR}
#-DDG_DEBUG -DDEBUG -DNOSCF  (MPI)

#> lib or other include file path
LIB_PATH_ROOT=../lib_gfortran
FFTW_DIR_ROOT=../fftw-3.3.9


LIBS = -L ${LIB_PATH_ROOT} -lfftw3 -lfmm3d -larpack_gfortran -lxcf90 -lxc -llapack -lrefblas
OMP_LIBS =-L ${LIB_PATH_ROOT} -lfftw3 -lfmm3d -larpack_gfortran -lxcf90 -lxc -lscalapack -llapack -lrefblas# -L/home/xlt/Software/mpiP-3.5/lib -lmpiP
OMP_FFT = -L ${FFTW_DIR_ROOT}/lib -lfftw3_mpi -lfftw3 
endif

#> ifort
ifeq ($(PLATFORM_ENV),'intel')
FC = ifort
LINKER = $(FC)
MPIFC = mpiifort
MPILINKER = $(MPIFC)

#FFLAGS =  -g -fpp -I${LIB_PATH_ROOT}/include -CA -CU -CB -traceback -check all -cpp -module ${WORKDIR}
FFLAGS = -O3 -fpp  -I${LIB_PATH_ROOT}/include  -module ${WORKDIR}
#OMP_FFLAGS = -g -CA -CB -CU -traceback -heap-arrays 1000 -fpp -DMPI -I$(FFTW_DIR_ROOT)/include -I${LIB_PATH_ROOT}/include -module ${OMP_WORKDIR}
OMP_FFLAGS = -O3 -g -traceback -heap-arrays 10 -align array32byte -align rec32byte -fpp -DMPI -I$(FFTW_DIR_ROOT)/include -I${LIB_PATH_ROOT}/include -module ${OMP_WORKDIR}
# -heap-arrays 1000
#> MACRO for FFLAGS(MPI)
#-DDG_DEBUG -DDEBUG -DNOSCF

#> lib or other include file path
LIB_PATH_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/ares_lib/lib
LIB_PATH=-L ${LIB_PATH_ROOT}
FFTW_DIR_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/ares_lib/lib/fftw
FFTW_DIR=-L ${FFTW_DIR_ROOT}/lib
LIB_MKL_ROOT=${MKLROOT}/lib/intel64
LIB_MKL_PATH=-L ${LIB_MKL_ROOT}


#LIBS = -lfftw3 -lfmm3d -larpack_ifort -lxcf90 -lxc
LIBS = -lfftw3 -larpack_ifort -lxcf90 -lxc
MKLLIBS= -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#OMP_LIBS = -lfmm3d  -larpack_ifort -lxcf90 -lxc# -L/home/xlt/Software/mpiP-3.5/lib -lmpiP
OMP_LIBS = -larpack_ifort -lxcf90 -lxc
OMP_FFT = -lfftw3_mpi -lfftw3 
OMP_MKLLIBS= -lmkl_blacs_intelmpi_lp64 -lmkl_lapack95_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5  -lpthread
endif

#----------------------------------------------------------------------------------
#> (2) PLATFORM INDEPENDENT VARIABLES.
#----------------------------------------------------------------------------------

GOAL = ares.x
OMP_GOAL = ares.mpi

MODULES = Constants.o  m_time_evaluate.o Parameters.o array_io.o MathSplines.o  Math.o  Fourier.o Struct_module.o Psp_module.o Read_module.o \
     Succeed_module.o \
	 Grid_module.o Ewald.o  Finite_module.o XC_functional.o \
     sgfft.o poisson_isf.o \
	 Begin_module.o Smearing_module.o Nonlocalpot_module.o \
	 IonLocalPotentialAssignment.o Isolate_module.o\
         Potential_module.o Lapack_module.o Matvec_module.o Arpack_module.o Chebyshev_fliter.o Energy_module.o \
	 Mixer_module.o Scf_module.o Out_module.o ForceStress_module.o\
     Write_module.o dyn_module.o brent.o cg_vasp.o cg_module.o \
     Relax_module.o End_module.o Bands_module.o

OMP_MODULES = Constants.o m_time_evaluate.o Parameters.o Smpi_math_module.o array_io.o MathSplines.o  Math.o  Fourier.o Struct_module.o Psp_module.o Read_module.o \
     Succeed_module.o \
     Scala.o \
	 Grid_module.o Ewald.o  Finite_module.o XC_functional.o \
     sgfft.o poisson_isf.o \
	 Begin_module.o Smearing_module.o Nonlocalpot_module.o \
	 IonLocalPotentialAssignment.o Isolate_module.o\
     Potential_module.o Lapack_module.o Matvec_module.o Arpack_module.o Chebyshev_fliter.o Energy_module.o \
	 Mixer_module.o Scf_module.o Out_module.o ForceStress_module.o \
     Write_module.o dyn_module.o brent.o cg_vasp.o cg_module.o \
     Relax_module.o End_module.o Bands_module.o

BASE = Main.o

#> containing object files for serial version
WORKDIR = ../bin
OMP_WORKDIR = ../ompbin


ALLS = $(MODULES) $(BASE)
OMP_ALLS = $(OMP_MODULES) $(BASE)

$(WORKDIR)/%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

# openmp
$(OMP_WORKDIR)/%.o : %.f90
	$(MPIFC) -c $(OMP_FFLAGS) $< -o $@

#> Absolute path of the object files (with working dir attached)
OBJS = $(patsubst %.o, $(WORKDIR)/%.o, $(ALLS))
OMP_OBJS = $(patsubst %.o, $(OMP_WORKDIR)/%.o, $(OMP_ALLS))

default :
	@ if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	@ chmod a+rwX $(WORKDIR)
	@ $(MAKE) $(GOAL)

mpi :
	@ if [ ! -d $(OMP_WORKDIR) ]; then mkdir $(OMP_WORKDIR); fi
	@ chmod a+rwX $(OMP_WORKDIR)
	@ $(MAKE) $(OMP_GOAL)
# Rules
#$(GOAL) : $(OBJS)
#	$(LINKER) $(FFLAGS) -o $(GOAL) $(OBJS)\
#        ../lapack-3.8.0/liblapack.a ../lapack-3.8.0/librefblas.a\
#     ${LIB_PATH} $(LIBS)
$(GOAL) : $(OBJS)
	$(LINKER) $(FFLAGS) -o $(GOAL) $(OBJS)\
     ${LIB_PATH} $(LIBS)\
     ${LIB_MKL_PATH} ${MKLLIBS}

#$(OMP_GOAL) : $(OMP_OBJS)
#	$(MPILINKER) $(OMP_FFLAGS) -o $(OMP_GOAL) $(OMP_OBJS)\
#     ${LIB_PATH} $(OMP_LIBS)\
#     ../scalapack-2.0.2/lib/libscalapack.a\
#        ../lapack-3.8.0/liblapack.a ../lapack-3.8.0/librefblas.a

$(OMP_GOAL) : $(OMP_OBJS)
	$(MPILINKER) $(OMP_FFLAGS) -o $(OMP_GOAL) $(OMP_OBJS)\
     ${FFTW_DIR} ${OMP_FFT}\
     ${LIB_PATH} $(OMP_LIBS)\
     ${LIB_MKL_PATH} ${OMP_MKLLIBS}

cleanall : clean cleanomp

all : default omp

.PHONY : cleanall all tags clean_mod

clean :
	@ if [ -d $(WORKDIR) ]; then $(RM) -r $(WORKDIR); fi
	@ if [ -f $(GOAL) ] ; then $(RM) $(GOAL); fi
cleanmpi :
	@ if [ -d $(OMP_WORKDIR) ]; then $(RM) -r $(OMP_WORKDIR); fi
	@ if [ -f $(OMP_GOAL) ]; then $(RM) $(OMP_GOAL); fi

clean_mod :
	@ if [ -d $(WORKDIR) ]; then $(RM) -r $(WORKDIR); fi

tags :
	@ ctags *.f90
.PHONY : relink
relink:
	@$(LINKER) $(FFLAGS) -o $(GOAL) $(OBJS) ${LIB_PATH} $(LIBS)  ${LIB_MKL_PATH} ${MKLLIBS}

tar_src :$(wildcard *90 *inc *.h Makefile )
	@tar -czvf ares_src.tar.gz $^

omprelink:
	$(MPILINKER) $(OMP_FFLAGS) -o $(OMP_GOAL) $(OMP_OBJS)\
     ${FFTW_DIR} ${OMP_FFT}\
        ../lapack-3.8.0 -lrefblas -llapack\
     ../scalapack-2.0.2/lib/libscalapack.a\
     ${LIB_PATH} $(OMP_LIBS)
