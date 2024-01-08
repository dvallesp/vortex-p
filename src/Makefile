NAME=vortex

FFTW=1 # default value: use FFTW; set to 0 if not installed 
$(info FFTW=$(FFTW))

ifeq ($(COMP),1)          
 FC=gfortran 
 FFLAGS=-O3 -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW)
 LIBS=-L/usr/local/lib/libcoretran.a -lcoretran
 INC=-I/usr/local/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/usr/local/lib64 -lfftw3f_omp -lfftw3f
 endif
endif

ifeq ($(COMP),999)
 FC=gfortran
 FFLAGS= -Og -fbounds-check -fbacktrace -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW)
 LIBS=-L/usr/local/lib/libcoretran.a -lcoretran
 INC=-I/usr/local/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/usr/local/lib64 -lfftw3f_omp -lfftw3f
 endif
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
