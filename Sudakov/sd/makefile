FFLAGS = -g -fno-automatic -ffixed-line-length-none  
FC=gfortran


OBJECTS =  main.o \
           vegas.o \
           alfas.o 

sdk:  $(OBJECTS)
	$(FC) -o ask $(FFLAGS) $(OBJECTS) 

f.o:  $(FC) -c -o $*.o  $(FFLAGS) $*.f
