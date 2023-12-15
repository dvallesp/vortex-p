NAME=vortex

ifeq ($(COMP),1)          
	FC=gfortran 
	FFLAGS=-O3 -fopenmp -mcmodel=medium 
	LIBS=-L/usr/local/lib/libcoretran.a -lcoretran
	INC=-I/usr/local/include/coretran/
endif

ifeq ($(COMP),999)
 FC=gfortran
 FFLAGS= -Og -fbounds-check -fbacktrace -fopenmp -mcmodel=medium
 LIBS=-L/usr/local/lib/libcoretran.a -lcoretran
 INC=-I/usr/local/include/coretran/
endif


##########################################################################

EXEC=$(NAME)

OBJ=read_gadget.o vortex.o

$(NAME): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $(EXEC) $(LIBS)

read_gadget.o: read_gadget.f 
	$(FC) $(FFLAGS) -c -o read_gadget.o read_gadget.f

vortex.o: vortex.f
	$(FC) $(FFLAGS) -static $(INC) -c -o vortex.o vortex.f

clean:
	 rm -f $(OBJ) $(EXEC)