.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

include ../../arch.make

ifneq (,$(findstring -DCDF,$(FPPFLAGS)))
# Only the save object has references to the netcdf library
INC += -I../../ncdf/src -I../../ncdf/fvar
endif

# The phony list
.PHONY: default lib deb clean

# We need to add the src-directory to the sources
FFLAGS += -I../../src

default: lib

# Make dependecy list
dep:
	sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90

# Grid objects
OBJS = m_mg_info.o m_time.o

LIB_TEST_NAME?=libmg_test.a
lib: $(LIB_TEST_NAME)
$(LIB_TEST_NAME): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB_TEST_NAME) $(OBJS)
	-$(RANLIB) $(LIB_TEST_NAME)

clean:
	rm -f $(OBJS) *.mod *.MOD $(LIB_TEST_NAME)

# DO NOT DELETE THIS LINE - used by make depend
