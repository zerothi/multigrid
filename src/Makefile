.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

include ../arch.make

# The phony list
.PHONY: default deb clean mg

ifeq ($(USE_CDF),1)
ifneq (,$(findstring -DMG__CDF,$(FPPFLAGS)))
FPPFLAGS += -DMG__CDF
endif
endif

ifneq (,$(findstring -DMG__CDF,$(FPPFLAGS)))
# Only the save object has references to the netcdf library
LIBS := ../ncdf/libncdf.a ../ncdf/fdict/libfdict.a $(LIBS)
INCLUDES += -I../ncdf -I../ncdf/fdict
endif

LIB_NAME?=libmg.a
default: $(LIB_NAME)

# Make dependecy list
dep:
	sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90

# Default objects
OBJS  = m_unit.o

# Grid objects
OBJS += t_mg.o t_mg_interp.o t_bc.o m_mg_save.o m_io.o m_mg_io.o

# Gauss-Seidel objects
OBJS += m_gs_CDS.o

# The application
mg: $(LIB_NAME) $(OBJS) mg.o
	$(FC) -o $@ mg.o $(OBJS) $(LIB_NAME) $(LIBS)

ifneq (,$(findstring -DCDF,$(FPPFLAGS)))
# Only the save object has references to the netcdf library
m_mg_save.o: INCLUDES+=-I../ncdf -I../ncdf/fdict
endif

# library
$(LIB_NAME): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB_NAME) $(OBJS)
	-$(RANLIB) $(LIB_NAME)

clean:
	rm -f $(LIB_NAME) $(OBJS) *.mod *.MOD

# DO NOT DELETE THIS LINE - used by make depend
m_gs_CDS.o: t_mg.o t_mg_interp.o
m_mg_io.o: m_io.o m_unit.o t_mg.o
m_mg_save.o: m_io.o t_mg.o
mg.o: m_gs_CDS.o m_io.o m_mg_io.o m_mg_save.o t_mg.o
t_mg.o: m_unit.o t_bc.o
t_mg_interp.o: t_mg.o
