.SUFFIXES: .f90 .o
#---------------------------------------------------------------------------------
#> (0)easy config
#---------------------------------------------------------------------------------
#> 'intel' or 'gcc'
#PLATFORM_ENV='gnu'
PLATFORM_ENV='intel'
#---------------------------------------------------------------------------------
#> (1)perpare for different platforms
#---------------------------------------------------------------------------------
#> containing object files for serial version
WORKDIR = ../bin
MPI_WORKDIR = ../mpibin
#> gfortran
ifeq ($(PLATFORM_ENV),'gcc')
FC = gfortran
LINKER = $(FC)
MPIFC = mpif90
MPILINKER = $(MPIFC)

FFLAGS = -O3 -cpp  -ffree-line-length-none  -I${LIB_PATH_ROOT}/include -J${WORKDIR}
FFLAGS_debug = -g -cpp  -ffree-line-length-none  -I${LIB_PATH_ROOT}/include -J${WORKDIR}
#FFLAGS = ${FFLAGS_debug}
MPI_FFLAGS = -g -cpp  -ffree-line-length-none  -Wtabs -fbounds-check -fcheck=all -DMPI -I${FFTW_DIR_ROOT}/include -I${LIB_PATH_ROOT}/include -J${MPI_WORKDIR}
#-fsanitize=address 
#MPI_FFLAGS = -g -O3 -cpp  -ffree-line-length-none  -fcheck=all -DMPI -I${LIB_PATH_ROOT}/include -I${FFTW_DIR_ROOT}/include -J${MPI_WORKDIR}
#-DDG_DEBUG -DDEBUG -DNOSCF  (MPI)

#> lib or other include file path
LIB_PATH_ROOT=../lib_gfortran
FFTW_DIR_ROOT=../fftw-3.3.9


LIBS = -L ${LIB_PATH_ROOT} -lfftw3 -lfmm3d -larpack_gfortran -lxcf90 -lxc -llapack -lrefblas
MPI_LIBS =-L ${LIB_PATH_ROOT} -lfftw3 -lfmm3d -larpack_gfortran -lxcf90 -lxc -lscalapack -llapack -lrefblas# -L/home/xlt/Software/mpiP-3.5/lib -lmpiP
MPI_FFT = -L ${FFTW_DIR_ROOT}/lib -lfftw3_mpi -lfftw3 
endif

#> ifort
ifeq ($(PLATFORM_ENV),'intel')
FC = ifort
LINKER = $(FC)
MPIFC = mpiifort
MPILINKER = $(MPIFC)

#> lib or other include file path
LIB_PATH_ROOT=/work/home/xinyu/soft/ARES/ares_lib/lib
LIB_PATH=-L ${LIB_PATH_ROOT}
LIB_MKL_ROOT=${MKLROOT}/lib/intel64
LIB_MKL_PATH=-L ${LIB_MKL_ROOT}
FFTW_DIR_ROOT=/work/home/xinyu/soft/ARES/ares_lib/lib/fftw
FFTW_DIR=-L ${FFTW_DIR_ROOT}/lib

#FFLAGS
FFLAGS = -O3 -fpp  -I${LIB_PATH_ROOT}/include  -module ${WORKDIR}
MPI_FFLAGS = -O3 -traceback -heap-arrays 1000 -align array32byte -align rec32byte -fpp -DMPI -I$(FFTW_DIR_ROOT)/include -I${LIB_PATH_ROOT}/include -module ${MPI_WORKDIR}
#Debug
#FFLAGS = -g -fpp  -I${LIB_PATH_ROOT}/include  -module ${WORKDIR}
#MPI_FFLAGS = -g -CU -CB -CA -traceback -heap-arrays 1000 -align array32byte -align rec32byte -fpp -DMPI -I$(FFTW_DIR_ROOT)/include -I${LIB_PATH_ROOT}/include -module ${MPI_WORKDIR}

#LIBS = -lfftw3 -lfmm3d -larpack_ifort -lxcf90 -lxc
LIBS = -lfftw3 -larpack_ifort -lxcf90 -lxc
MKLLIBS= -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#MPI_LIBS = -lfmm3d  -larpack_ifort -lxcf90 -lxc# -L/home/xlt/Software/mpiP-3.5/lib -lmpiP
MPI_LIBS = -larpack_ifort -lxcf90 -lxc
MPI_FFT = -lfftw3_mpi -lfftw3 
MPI_MKLLIBS= -lmkl_blacs_intelmpi_lp64 -lmkl_lapack95_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5  -lpthread
endif

#----------------------------------------------------------------------------------
#> (2) PLATFORM INDEPENDENT VARIABLES.
#----------------------------------------------------------------------------------

GOAL = rsdft.x
MPI_GOAL = ares.pbc

MODULES = Constants.o  m_time_evaluate.o Parameters.o array_io.o MathSplines.o  Math.o  \
	Fourier.o Struct_module.o Psp_module.o Read_module.o \
     Succeed_module.o \
	 Grid_module.o Ewald.o  Finite_module.o XC_functional.o \
     sgfft.o poisson_isf.o \
	 Begin_module.o Smearing_module.o Nonlocalpot_module.o \
	 IonLocalPotentialAssignment.o Isolate_module.o\
         Potential_module.o Lapack_module.o Matvec_module.o Arpack_module.o Chebyshev_fliter.o Energy_module.o \
	 Mixer_module.o Scf_module.o Out_module.o ForceStress_module.o\
     Write_module.o dyn_module.o brent.o cg_vasp.o cg_module.o \
     Relax_module.o End_module.o Bands_module.o

MPI_MODULES = Constants.o MPI_time_evaluate.o Parameters.o MPI_array_io.o MathSplines.o  Math.o \
	Smpi_math_module.o Fourier.o Struct_module.o Read_module.o \
	Succeed_module.o ScaLapack_module.o Ewald.o  Finite_module.o XC_functional.o Nonlocalpot_module.o \
	Begin_module.o Smearing_module.o \
	Potential_module.o Lapack_module.o Matvec_module.o Arpack_module.o Chebyshev_fliter.o Energy_module.o Mixer_module.o Scf_module.o 

BASE = Main.o

#> containing object files for serial version
WORKDIR = ../bin
MPI_WORKDIR = ../mpibin


ALLS = $(MODULES) $(BASE)
MPI_ALLS = $(MPI_MODULES) $(BASE)

$(WORKDIR)/%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

# mpi
$(MPI_WORKDIR)/%.o : %.f90
	$(MPIFC) -c $(MPI_FFLAGS) $< -o $@

#> Absolute path of the object files (with working dir attached)
OBJS = $(patsubst %.o, $(WORKDIR)/%.o, $(ALLS))
MPI_OBJS = $(patsubst %.o, $(MPI_WORKDIR)/%.o, $(MPI_ALLS))

default :
	@ if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	@ chmod a+rwX $(WORKDIR)
	@ $(MAKE) $(GOAL)

# mpi :
# 	@ if [ ! -d $(MPI_WORKDIR) ]; then mkdir $(MPI_WORKDIR); fi
# 	@ chmod a+rwX $(MPI_WORKDIR)
# 	@ $(MAKE) $(MPI_GOAL)
# Rules
#$(GOAL) : $(OBJS)
#	$(LINKER) $(FFLAGS) -o $(GOAL) $(OBJS)\
#        ../lapack-3.8.0/liblapack.a ../lapack-3.8.0/librefblas.a\
#     ${LIB_PATH} $(LIBS)
$(GOAL) : $(OBJS)
	$(LINKER) $(FFLAGS) -o $(GOAL) $(OBJS)\
     ${LIB_PATH} $(LIBS)\
     ${LIB_MKL_PATH} ${MKLLIBS}


$(MPI_GOAL) : $(MPI_OBJS)
	$(MPILINKER) $(MPI_FFLAGS) -o $(MPI_GOAL) $(MPI_OBJS)\
     ${FFTW_DIR} ${MPI_FFT}\
     ${LIB_PATH} $(MPI_LIBS)\
     ${LIB_MKL_PATH} ${MPI_MKLLIBS}

cleanall : clean cleanmpi

all : default mpi

.PHONY : cleanall all tags clean_mod

clean :
	@ if [ -d $(WORKDIR) ]; then $(RM) -r $(WORKDIR); fi
	@ if [ -f $(GOAL) ] ; then $(RM) $(GOAL); fi
cleanmpi :
	@ if [ -d $(MPI_WORKDIR) ]; then $(RM) -r $(MPI_WORKDIR); fi
	@ if [ -f $(MPI_GOAL) ]; then $(RM) $(MPI_GOAL); fi

clean_mod :
	@ if [ -d $(WORKDIR) ]; then $(RM) -r $(WORKDIR); fi

tags :
	@ ctags *.f90
.PHONY : relink
relink:
	@$(LINKER) $(FFLAGS) -o $(GOAL) $(OBJS) ${LIB_PATH} $(LIBS)  ${LIB_MKL_PATH} ${MKLLIBS}

tar_src :$(wildcard *90 *inc *.h Makefile )
	@tar -czvf ares_src.tar.gz $^

mpirelink:
	$(MPILINKER) $(MPI_FFLAGS) -o $(MPI_GOAL) $(MPI_OBJS)\
     ${FFTW_DIR} ${MPI_FFT}\
        ../lapack-3.8.0 -lrefblas -llapack\
     ../scalapack-2.0.2/lib/libscalapack.a\
     ${LIB_PATH} $(MPI_LIBS)
