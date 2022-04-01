Options_for_Debugging=-o0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interface -warn interface

Flag=$(Options_for_Debugging)

FC=ifort

targets=3_1


build:
	$(FC) -c  $(targets).f90  $(Flag)
	$(FC) -o  $(targets) $(targets).o 

