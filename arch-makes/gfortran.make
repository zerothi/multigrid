FC=gfortran
FC_SERIAL=gfortran
AR=ar
RANLIB=ranlib

FFLAGS=-O3 -m64 -fPIC -fno-second-underscore -ftree-vectorize -fexpensive-optimizations \
	-finline-functions-called-once -funroll-loops -fvariable-expansion-in-unroller \
	-ftree-loop-optimize -frename-registers -fprefetch-loop-arrays -finline-small-functions \
	-fipa-pure-const -foptimize-sibling-calls -fipa-cp

FFLAGS += -fopenmp

#FFLAGS = -g -O0 -Warray-bounds

#FFLAGS += -Wunused

ARFLAGS = cru

C_V=gnu-4.9.2
Z_PATH=/opt/zlib/1.2.8/$(C_V)
NCDF_PATH=/opt/netcdf-serial/4.3.3/$(C_V)
HDF5_PATH=/opt/hdf5-serial/1.8.14/$(C_V)

INC = -I$(NCDF_PATH)/include

FPPFLAGS = -DCDF

LIBS = -fopenmp

ifneq (,$(findstring -DCDF,$(FPPFLAGS)))
LIBS += -L$(NCDF_PATH)/lib -Wl,-rpath=$(NCDF_PATH)/lib
LIBS += -L$(HDF5_PATH)/lib -Wl,-rpath=$(HDF5_PATH)/lib
LIBS += -L$(Z_PATH)/lib -Wl,-rpath=$(Z_PATH)/lib
LIBS += -lnetcdff -lnetcdf -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz
endif

.F90.o:
	$(FC) -c $(INC) $(FFLAGS) $(FPPFLAGS) $< 
.f90.o:
	$(FC) -c $(INC) $(FFLAGS) $<

