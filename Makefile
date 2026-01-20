CMP     = intel
OPENMP  = 1

ifeq ($(OPENMP),1)
  DEFS_ADD += -DOPENMP
endif

exeName =  cd
dir = ./src/
inc = 
lib = 

ifeq ($(CMP),intel)
  FortC = mpiifx
  CFLAG = -fpp -O3 -mcmodel=large -fp-model fast=2 -g#-heap-arrays 20#-ipo -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2
else ifeq ($(CMP),gcc)
  FortC = mpif90
  CFLAG = -cpp -O3 -Warray-bounds -fbacktrace -g -ffree-line-length-none -mcmodel=large -funroll-loops -floop-optimize
endif

SrcT   := m_parameters.f90 healpix_types.f90 bit_manipulation.f90 m_Healpix.f90 m_Raytheia.f90 m_readdensity.f90 m_calc_columndens.f90 m_outputs.f90 main.f90
src:= $(addprefix $(dir), ${SrcT})

all: $(exeName)
src_obj  = $(src:%.f90=%.o)
$(exeName):$(src_obj)
	$(FortC) $(CFLAG) -fopenmp -o $@ $(src_obj) $(lib)
$(src_obj):%.o :%.f90
	$(FortC) $(CFLAG) -fopenmp $(inc) $(DEFS_ADD) -c $<
	@ mv $(@F) ${dir}

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(dir)*.o

