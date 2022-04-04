ARESSRC=Arpack_module.f90 Bands_module.f90 Begin_module.f90 \
		Chebyshev_fliter.f90 Constants.f90 End_module.f90 Energy_module.f90 \
		Ewald.f90 Finite_module.f90 ForceStress_module.f90 Fourier.f90 \
		Lapack_module.f90 Math.f90 MathSplines.f90 Matvec_module.f90 \
		Mixer_module.f90 MPI_array_io.f90 MPI_time_evaluate.f90 \
		Nonlocalpot_module.f90 Parameters.f90 Potential_module.f90 \
		Read_module.f90 ScaLapack_module.f90 Scf_module.f90 \
		Smearing_module.f90 Smpi_math_module.f90 Struct_module.f90 \
		Succeed_module.f90 Write_module.f90 XC_functional.f90

MPI_MODULES = Constants.o MPI_time_evaluate.o Parameters.o MPI_array_io.o \
	MathSplines.o Math.o \
	Smpi_math_module.o Fourier.o Struct_module.o Read_module.o \
	Succeed_module.o ScaLapack_module.o Ewald.o  Finite_module.o \
	XC_functional.o Nonlocalpot_module.o \
	Begin_module.o Smearing_module.o \
	Potential_module.o Lapack_module.o Matvec_module.o \
	Arpack_module.o Chebyshev_fliter.o Energy_module.o Mixer_module.o Scf_module.o 



F90WRAP_FILES=f90wrap_Arpack_module.f90 f90wrap_Finite_module.f90 f90wrap_MPI_array_io.f90 \
				f90wrap_Smearing_module.f90 f90wrap_Bands_module.f90 f90wrap_ForceStress_module.f90 \
				f90wrap_MPI_time_evaluate.f90 f90wrap_Smpi_math_module.f90 f90wrap_Begin_module.f90 \
				f90wrap_Fourier.f90 f90wrap_Nonlocalpot_module.f90 f90wrap_Struct_module.f90 \
				f90wrap_Chebyshev_fliter.f90 f90wrap_Lapack_module.f90 f90wrap_Parameters.f90 \
				f90wrap_Succeed_module.f90 f90wrap_Constants.f90 f90wrap_Math.f90 \
				f90wrap_Potential_module.f90 f90wrap_Write_module.f90 f90wrap_End_module.f90 \
				f90wrap_MathSplines.f90 f90wrap_Read_module.f90 f90wrap_XC_functional.f90 \
				f90wrap_Energy_module.f90 f90wrap_Matvec_module.f90 f90wrap_ScaLapack_module.f90 \
				f90wrap_Ewald.f90 f90wrap_Mixer_module.f90 f90wrap_Scf_module.f90



FC = ifort
LINKER = $(FC)
MPIFC = mpiifort
MPILINKER = $(MPIFC)
LIB_PATH_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/ares_lib/lib
LIB_PATH=-L${LIB_PATH_ROOT}
LIB_MKL_ROOT=${MKLROOT}/lib/intel64
LIB_MKL_PATH=-L${LIB_MKL_ROOT}
FFTW_DIR_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/ares_lib/lib/fftw
FFTW_DIR=-L${FFTW_DIR_ROOT}/lib

FFLAGS = -O3 -fpp  -I${LIB_PATH_ROOT}/include  -module .
MPI_FFLAGS = -O3 -traceback -heap-arrays 1000 -align array32byte -align rec32byte -fpp -DMPI -I$(FFTW_DIR_ROOT)/include -I${LIB_PATH_ROOT}/include -module .

LIBS = -lfftw3 -larpack_ifort -lxcf90 -lxc
MKLLIBS= -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
MPI_LIBS = -larpack_ifort -lxcf90 -lxc
MPI_FFT = -lfftw3_mpi -lfftw3 
MPI_MKLLIBS= -lmkl_blacs_intelmpi_lp64 -lmkl_lapack95_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5  -lpthread

ifneq ($(findstring -std, ${CFLAGS}), -std)
   CFLAGS := $(CFLAGS) -std=c99
endif


mpiWrap:
	@$(MAKE) f90WrapFile
	@$(MAKE) buildObj
	CFLAGS='$(CFLAGS)' \
	f2py --fcompiler=intelem \
		-c --f90exec=$(MPIFC) --f77exec=$(MPIFC) --f90flags='${MPI_FFLAGS}' \
		-m _arespy ${F90WRAP_FILES} $(MPI_MODULES) ${FFTW_DIR} ${MPI_FFT} \
		${LIB_PATH} $(MPI_LIBS) ${LIB_MKL_PATH} ${MPI_MKLLIBS}

f90WrapFile:
	f90wrap -m arespy $(ARESSRC) -k kindMap.json

buildObj:$(MPI_MODULES)
	echo $(MPI_MODULES)

$(MPI_MODULES): %.o: %.f90
	$(MPIFC) -c $(MPI_FFLAGS) $< -o $@

test :
	echo $(ARESSRC)