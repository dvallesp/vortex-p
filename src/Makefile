NAME=vortex

# whether to use FFTW; default value: use FFTW (1); set to 0 if not installed 
FFTW := $(if $(FFTW),$(FFTW),1)
$(info FFTW=$(FFTW))

# whether to use the turbulent filter; default value: no (0); set to 1 if filtering is desired
FILTER := $(if $(FILTER),$(FILTER),0)
$(info FILTER=$(FILTER))


#### This is an example for a local machine ####
ifeq ($(COMP),1)          
 FC=gfortran 
 FFLAGS=-O3 -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW) -Duse_filter=$(FILTER)
 LIBS=/usr/local/lib/libcoretran.a -lcoretran
 INC=-I/usr/local/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/usr/local/lib64 -lfftw3f_omp -lfftw3f
 endif
endif

#### This is for the LluisVives cluster at Universitat de València, where I locally installed coretran ####
ifeq ($(COMP),2)
 FC=gfortran
 FFLAGS=-O3 -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW) -Duse_filter=$(FILTER)
 LIBS=-Wl,-rpath,$(HOME)/coretran_shared_gnu/coretran/lib $(HOME)/coretran_shared_gnu/coretran/lib/libcoretran.so
 INC=-I$(HOME)/coretran_shared_gnu/coretran/library/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/storage/apps/FFTW/3.3.8/INTEL/IMPI/lib -lfftw3f_omp -lfftw3f
 endif
endif

#### Debugging option in my local machine ####
ifeq ($(COMP),999)
 FC=gfortran
 FFLAGS= -Og -fbounds-check -fbacktrace -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW) -Duse_filter=$(FILTER)
 LIBS=/usr/local/lib/libcoretran.a -lcoretran
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

info:
	@echo ""
	@echo "***********************************************************"
	@echo "***                       vortex-p                      ***"
	@echo "***********************************************************"
	@echo "* David Vallés-Pérez et al., Universitat de València, 2024*"
	@echo "***********************************************************"
	@echo "*** FLAGS for compiling ***"
	@echo "- COMP (mandatory): chooses the compiler options; provided "
	@echo "                    are three examples, but you might need " 
	@echo "                    to modify them for your own machine(s)."
	@echo "   In particular, you might need to change the paths to the "
	@echo "     coretran library and to the FFTW library (if using it)"
	@echo ""
	@echo "- FFTW (optional, default is 1): set to 0 if you do not have "
	@echo "                                 FFTW"
	@echo ""
	@echo "- FILTER (optional, default is 0): set to 1 if you want to "
	@echo "                                   use the turbulent filter"
	@echo ""
