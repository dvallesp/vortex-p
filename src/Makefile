NAME=vortex

ifeq ($(COMP),1)          
	FC=gfortran 
	FFLAGS=-O3 -fopenmp -mcmodel=medium 
	LIBS=libfortran_kdtree.a
endif

ifeq ($(COMP),999)
 FC=gfortran
 FFLAGS= -Og -fbounds-check -fbacktrace -fopenmp -mcmodel=medium
 LIBS=libfortran_kdtree.a
endif


##########################################################################

EXEC=$(NAME)

OBJ=read_gadget.o vortex.o


$(NAME): $(OBJ)

	$(FC) $(FFLAGS) $(OBJ) -o $(EXEC) $(LIBS) 

clean:
	 rm -f $(OBJ) $(EXEC)
