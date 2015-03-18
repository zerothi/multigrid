FC=ifort
FC_SERIAL=ifort
AR=xiar
RANLIB=ranlib

# -warn unused
# For debugging:
FFLAGS= -O0 -g -check bounds -traceback -fp-model strict -I$(INTEL_PATH)/mkl/include -I$(INTEL_PATH)/mkl/include/intel64/lp64 -mkl=sequential $(INC_PATH) 

#FFLAGS=-O1 -fp-model strict -I$(INTEL_PATH)/mkl/include -I$(INTEL_PATH)/mkl/include/intel64/lp64 -mkl=sequential $(INC_PATH)

FFLAGS += -openmp

LIBS = -openmp

ARFLAGS = cru

.F90.o:
	$(FC) -c $(FFLAGS) $(FPPFLAGS) $< 
.f90.o:
	$(FC) -c $(FFLAGS) $<

