Options_for_Debugging=-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interface -warn interface

WORKDIR=.

FFLAGS=$(Options_for_Debugging)

FC=ifort

targets=6_1

build:$(targets).o
	$(FC) $(FFLAGS) -o  $(targets) $<

$(WORKDIR)/%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@
