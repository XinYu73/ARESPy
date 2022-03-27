FC = ifort

FFLAGS = -fPIC

SRCS_F = $(wildcard matMult_*.f90)
OBJS   = $(SRCS_F:.f90=.o)

all: $(OBJS)
	ar crs libmatMultiplication.a $(OBJS)

clean:
	rm -f *.o  *.mod  *.so  *.a


.SUFFIXES:
.SUFFIXES: .o .f .f90 .h

.F90.o: $(FSOURCE)
	$(FC) -c $(FFLAGS) $?
