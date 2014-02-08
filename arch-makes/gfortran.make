FC=gfortran
FC_SERIAL=gfortran
AR=ar
RANLIB=ranlib

FFLAGS=-O3 -m64 -fPIC -fno-second-underscore -ftree-vectorize -fexpensive-optimizations \
	-finline-functions-called-once -funroll-loops -fvariable-expansion-in-unroller \
	-ftree-loop-optimize -frename-registers -fprefetch-loop-arrays -finline-small-functions \
	-fipa-pure-const -foptimize-sibling-calls -fipa-cp

FFLAGS += -fopenmp

#FFLAGS = -g -O0

#FFLAGS += -Wunused

ARFLAGS = cru

.F90.o:
	$(FC) -c $(FFLAGS) $(FPPFLAGS) $< 
.f90.o:
	$(FC) -c $(FFLAGS) $<

