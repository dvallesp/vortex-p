NAME=vortex

# whether to use FFTW; default value: use FFTW (1); set to 0 if not installed 
FFTW := $(if $(FFTW),$(FFTW),1)
$(info FFTW=$(FFTW))

# whether to use the turbulent filter; default value: no (0); set to 1 if filtering is desired
FILTER := $(if $(FILTER),$(FILTER),0)
$(info FILTER=$(FILTER))

# weighting scheme: number-weighted (0), mass-weighted (1), or volume-weighted (2)
WEIGHT := $(if $(WEIGHT),$(WEIGHT),0)
$(info WEIGHT=$(WEIGHT))

# kernel family (0: cubic spline; 1: Wendland C4; 2: Wendland C6)
KERNEL := $(if $(KERNEL),$(KERNEL),0)
$(info KERNEL=$(KERNEL))

# output preferences (for efficient compilation)
OUTPUT_GRID := $(if $(OUTPUT_GRID),$(OUTPUT_GRID),1)
$(info OUTPUT_GRID=$(OUTPUT_GRID))

OUTPUT_PARTICLES := $(if $(OUTPUT_PARTICLES),$(OUTPUT_PARTICLES),1)
$(info OUTPUT_PARTICLES=$(OUTPUT_PARTICLES))

OUTPUT_FILTER:=$(if $(OUTPUT_FILTER),$(OUTPUT_FILTER),1)
ifeq ($(FILTER),0)
  OUTPUT_FILTER=0
endif
$(info OUTPUT_FILTER=$(OUTPUT_FILTER))

#### This is an example for a local machine ####
ifeq ($(COMP),1)          
 FC=gfortran 
 FFLAGS=-O3 -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW) -Duse_filter=$(FILTER) -Dweight_scheme=$(WEIGHT) -Doutput_grid=$(OUTPUT_GRID) -Doutput_particles=$(OUTPUT_PARTICLES) -Doutput_filter=$(OUTPUT_FILTER) -Dikernel=$(KERNEL)
 LIBS=/usr/local/lib/libcoretran.a -lcoretran
 INC=-I/usr/local/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/usr/local/lib64 -lfftw3f_omp -lfftw3f
 endif
endif

#### This is for the LluisVives cluster at Universitat de València, where I locally installed coretran ####
ifeq ($(COMP),2)
 FC=gfortran
 FFLAGS=-O3 -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW) -Duse_filter=$(FILTER) -Dweight_scheme=$(WEIGHT) -Doutput_grid=$(OUTPUT_GRID) -Doutput_particles=$(OUTPUT_PARTICLES) -Doutput_filter=$(OUTPUT_FILTER) -Dikernel=$(KERNEL)
 LIBS=-Wl,-rpath,$(HOME)/coretran_shared_gnu/coretran/lib $(HOME)/coretran_shared_gnu/coretran/lib/libcoretran.so
 INC=-I$(HOME)/coretran_shared_gnu/coretran/library/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/storage/apps/FFTW/3.3.8/INTEL/IMPI/lib -lfftw3f_omp -lfftw3f
 endif
endif

#################################################################################
#### This is for the installation of coretran suggested in the documentation ####
ifeq ($(COMP),3)
 FC=gfortran
 FFLAGS=-O3 -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW) -Duse_filter=$(FILTER) -Dweight_scheme=$(WEIGHT) -Doutput_grid=$(OUTPUT_GRID) -Doutput_particles=$(OUTPUT_PARTICLES) -Doutput_filter=$(OUTPUT_FILTER) -Dikernel=$(KERNEL)
 LIBS=-Wl,-rpath,$(HOME)/coretran/lib $(HOME)/coretran/lib/libcoretran.so
 INC=-I$(HOME)/coretran/library/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/usr/local/lib64 -lfftw3f_omp -lfftw3f
 endif
endif
#################################################################################

#### Debugging option in my local machine ####
ifeq ($(COMP),999)
 FC=gfortran
 FFLAGS= -Og -fbounds-check -fbacktrace -fopenmp -mcmodel=medium -cpp -Duse_fftw=$(FFTW) -Duse_filter=$(FILTER) -Dweight_scheme=$(WEIGHT) -Doutput_grid=$(OUTPUT_GRID) -Doutput_particles=$(OUTPUT_PARTICLES) -Doutput_filter=$(OUTPUT_FILTER) -Dikernel=$(KERNEL)
 LIBS=/usr/local/lib/libcoretran.a -lcoretran
 INC=-I/usr/local/include/coretran/
 ifeq ($(FFTW),1)
  LIBS +=-L/usr/local/lib64 -lfftw3f_omp -lfftw3f
 endif
endif


##########################################################################

EXEC=$(NAME)

OBJ=read_gadget.o particle_data.o vortex.o

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
	@echo "- FFTW (optional, default: 1): set to 0 if you do not have "
	@echo "                    FFTW"
	@echo ""
	@echo "- FILTER (optional, default: 0): set to 1 if you want to use "
	@echo "                    the turbulent filter"
	@echo ""
	@echo "- WEIGHT (optional, default: 0): 0: number-weighted"
	@echo "                                 1: mass-weighted"
	@echo "                                 2: volume-weighted"	
	@echo ""
	@echo "- KERNEL (optional, default: 0): 0: cubic spline"
	@echo "                                 1: Wendland C4"
	@echo "                                 2: Wendland C6"
	@echo ""
	@echo "- OUTPUT_GRID (optional, default: 1): set to 0 if you do "
	@echo "                    not want to output gridded data"
	@echo ""
	@echo "- OUTPUT_PARTICLES (optional, default: 1): set to 0 if you "
	@echo "                    do not want to output particle data"
	@echo ""
	@echo "- OUTPUT_FILTER (optional, default: 1): set to 0 if you do "
	@echo "                    not want to output data about the "
	@echo "                    turbulent filter"
	@echo ""
