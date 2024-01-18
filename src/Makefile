NAME=vortex

# FFTW default value: use FFTW; set to 0 if not installed 
FFTW := $(if $(FFTW),$(FFTW),1)
$(info FFTW=$(FFTW))

#### This is an example for a local machine ####
ifeq ($(COMP),1)          
 FC=gfortran 
 FFLAGS=-O3 -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW)
 LIBS=/usr/local/lib/libcoretran.a -lcoretran
 INC=-I/usr/local/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/usr/local/lib64 -lfftw3f_omp -lfftw3f
 endif
endif

#### This is for the LluisVives cluster at Universitat de València, where I locally installed coretran ####
ifeq ($(COMP),2)
 FC=ifort
 FFLAGS=-O3 -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW)
 LIBS=-Wl,-rpath,$(HOME)/coretran_shared/coretran/lib $(HOME)/coretran_shared/coretran/lib/libcoretran.so
 INC=-I$(HOME)/coretran_shared/coretran/library/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/storage/apps/FFTW/3.3.8/INTEL/IMPI/lib -lfftw3f_omp -lfftw3f
 endif
endif

#### Debugging option in my local machine ####
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
	$(FC) $(FFLAGS) $(INC) -c -o vortex.o vortex.f

clean:
	 rm -f $(OBJ) $(EXEC)
