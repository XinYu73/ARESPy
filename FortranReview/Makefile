Options_for_Debugging=-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interface -warn interface

WORKDIR=.

FFLAGS=$(Options_for_Debugging)

FC=ifort

targets=7_2
objs=7_21.o 7_22.o
$(targets):$(objs)
	echo $^
	$(FC) $(FFLAGS) -o  $(targets) $(objs)

$(objs): %.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

.PHONY:clean
clean:
	rm -f $(targets) $(objs)
