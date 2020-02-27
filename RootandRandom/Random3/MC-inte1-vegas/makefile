FFLAGS = -g -m64 -fno-automatic -ffixed-line-length-none  
FC=gfortran

OBJECTS =  vegas.o \
           main.o \

lq:  $(OBJECTS)
	$(FC) $(FFLAGS) -o lq $(OBJECTS) 

f.o:  $(FC) -c -o $*.o  $(FFLAGS) $*.f
