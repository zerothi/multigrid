.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

include ../arch.make

# The phony list
.PHONY: default deb clean

LIB_NAME?=libmg.a
default: $(LIB_NAME)

# Make dependecy list
dep:
	sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90

# Grid objects
OBJS = t_mg.o m_mg_save.o

# Gauss-Seidel objects
OBJS += m_gs_CDS.o


# library
$(LIB_NAME): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB_NAME) $(OBJS)
	-$(RANLIB) $(LIB_NAME)

clean:
	rm -f $(LIB_NAME) $(OBJS) *.mod *.MOD

# DO NOT DELETE THIS LINE - used by make depend
m_gs_CDS.o: t_mg.o
m_gs_FOS.o: t_mg.o
m_gs.o: m_gs_FOS.o
