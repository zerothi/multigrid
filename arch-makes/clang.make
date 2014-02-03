FC=clang
FC_SERIAL=clang
AR=ar
RANLIB=ranlib

FFLAGS=-O3 -m64 -fPIC -funroll-loops -freroll-loops \
	$(INC_PATH)

#FFLAGS += -Wunused

ARFLAGS = cru

.F90.o:
	$(FC) -c $(FFLAGS) $(FPPFLAGS) $< 
.f90.o:
	$(FC) -c $(FFLAGS) $<

