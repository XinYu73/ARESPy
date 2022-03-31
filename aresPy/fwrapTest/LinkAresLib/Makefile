MyLIB_PATH_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/fwrapTest/f2pyStaticLib/testlib
MyLIBS=-L${MyLIB_PATH_ROOT} -lmatMultiplication

AresLib1=-L/work/home/xinyu/soft/ARES/ares_lib/lib  -larpack_ifort -lxcf90 -lxc
AresLib2=-L/work/home/xinyu/soft/ARES/ares_lib/lib/fftw/lib -lfftw3
AresLibInclude1=-I/work/home/xinyu/soft/ARES/ares_lib/lib/include
AresLibInclude2=-I/work/home/xinyu/soft/ARES/ares_lib/lib/fftw/include

LAPACKLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/lapack-3.10.0 -llapack -lrefblas  -ltmglib
BLASLIB=-L/work/home/xinyu/workplace/PhdProgram/LAPACK_BLAS/BLAS-3.10.0 -lblas


FC=ifort

.PHONY: buildLib Wrap
buildLib: $(MyLIB_PATH_ROOT)/*.o
	ar rc $(MyLIB_PATH_ROOT)/libmatMultiplication.a $(LIB_PATH_ROOT)/*.o

$(LIB_PATH_ROOT)/%.o:%.f90 
	$(FC) -c -fPIC $<

Wrap:
	@#f90wrap -m my my.f90
	@#f2py-f90wrap --f90exec=$(FC) --f77exec=$(FC) -c my.f90 f90wrap_toplevel.f90 -m my $(MyLIBS) $(LAPACKLIB) $(BLASLIB)

	
	@#f2py --f90exec=$(FC) --f77exec=$(FC) -c my.f90 -m my $(MyLIBS) $(LAPACKLIB) $(BLASLIB)

	f2py --f90exec=$(FC) --f77exec=$(FC) -c my.f90 -m my 

build:
	$(FC) -c  test.f90
	#$(FC) -o test test.o $(MyLIBS) $(AresLibInclude1) $(AresLibInclude2) $(AresLib1) $(AresLib2) 
	$(FC) -o test test.o -mkl