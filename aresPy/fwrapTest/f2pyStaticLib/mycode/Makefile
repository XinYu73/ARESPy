LIB_PATH_ROOT=/work/home/xinyu/workplace/PhdProgram/aresPy/fwrapTest/f2pyStaticLib/testlib
FC=ifort

LIBS=-L${LIB_PATH_ROOT} -lmatMultiplication

.PHONY: buildLib Wrap
buildLib: $(LIB_PATH_ROOT)/*.o
	ar rc $(LIB_PATH_ROOT)/libmatMultiplication.a $(LIB_PATH_ROOT)/*.o

$(LIB_PATH_ROOT)/%.o:%.f90 
	FC -c -fPIC $<

Wrap:
	f2py --f90exec=$(FC) --f77exec=$(FC) $(LIBS) -c matrixMult.f90 -m matrixMult