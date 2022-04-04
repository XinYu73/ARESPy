MyLIB_PATH_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/fwrapTest/f2pyStaticLib/testlib
MyLIBS=-L${MyLIB_PATH_ROOT} -lmatMultiplication

AresLib1=-L/work/home/xinyu/soft/ARES/ares_lib/lib  -larpack_ifort -lxcf90 -lxc
AresLib2=-L/work/home/xinyu/soft/ARES/ares_lib/lib/fftw/lib -lfftw3
AresLibInclude1=-I/work/home/xinyu/soft/ARES/ares_lib/lib/include
AresLibInclude2=-I/work/home/xinyu/soft/ARES/ares_lib/lib/fftw/include

LAPACKLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/lapack-3.10.0 -llapack -lrefblas  -ltmglib
BLASLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/BLAS-3.10.0 -lblas

Options_for_Debugging=-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interface -warn interface

WORKDIR=.

FFLAGS=$(Options_for_Debugging)

FC=ifort

.PHONY: buildLib Wrap
buildLib: $(MyLIB_PATH_ROOT)/*.o
	ar rc $(MyLIB_PATH_ROOT)/libmatMultiplication.a $(LIB_PATH_ROOT)/*.o

$(LIB_PATH_ROOT)/%.o:%.f90 
	$(FC) -c -fPIC $<

Wrap:
	@f2py Integrate.f90 -m Integrate -h Integrate.pyf --overwrite-signature
	@python -m numpy.f2py --fcompiler=intelem --f90exec=$(FC) --f77exec=$(FC) -c Integrate.pyf Integrate.f90 -m Integrate $(MyLIBS) $(LAPACKLIB) $(BLASLIB)

testobjs=AlgebraMethod.o Integrate.o test.o
test:$(testobjs)
	$(FC) -o test $^ $(LAPACKLIB) $(BLASLIB)
test1objs=AlgebraMethod.o Integrate.o test1.o 
test1:$(test1objs)
	$(FC) -o test1 $^ $(LAPACKLIB) $(BLASLIB)

$(test1objs): %.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(testobjs): %.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@