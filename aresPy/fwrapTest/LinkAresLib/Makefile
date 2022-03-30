MyLIB_PATH_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/fwrapTest/f2pyStaticLib/testlib
MyLIBS=-L${MyLIB_PATH_ROOT} -lmatMultiplication

AresLib1=-L/work/home/xinyu/soft/ARES/ares_lib/lib  -larpack_ifort -lxcf90 -lxc
AresLib2=-L/work/home/xinyu/soft/ARES/ares_lib/lib/fftw/lib -lfftw3
AresLibInclude1=-I/work/home/xinyu/soft/ARES/ares_lib/lib/include
AresLibInclude2=-I/work/home/xinyu/soft/ARES/ares_lib/lib/fftw/include

FC=ifort

.PHONY: buildLib Wrap
buildLib: $(MyLIB_PATH_ROOT)/*.o
	ar rc $(MyLIB_PATH_ROOT)/libmatMultiplication.a $(LIB_PATH_ROOT)/*.o

$(LIB_PATH_ROOT)/%.o:%.f90 
	$(FC) -c -fPIC $<

Wrap:
	f2py --f90exec=$(FC) --f77exec=$(FC) --opt=-O3 $(MyLIBS) $(AresLibInclude1) $(AresLibInclude2) $(AresLib1) $(AresLib2) -c my.f90 -m my 

build:
	$(FC) -c  test.f90
	$(FC) -o test test.o $(MyLIBS) $(AresLibInclude1) $(AresLibInclude2) $(AresLib1) $(AresLib2) 