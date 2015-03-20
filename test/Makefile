.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

# arch.make for correcting input
include ../arch.make

ifneq (,$(findstring -DCDF,$(FPPFLAGS)))
# Only the save object has references to the netcdf library
LIBS := ../ncdf/src/libncdf.a ../ncdf/fvar/libvardict.a $(LIBS)
INC += -I../ncdf/src -I../ncdf/fvar
endif

# The phony list
.PHONY: default lib deb clean all $(TESTS)

# We need to add the src-directory to the sources
FFLAGS += -I../src
INC += -Ilib

# Define the libraries
LIB_NAME?=libmg.a
LIB_NAME:=../src/$(LIB_NAME)
LIB_TEST_NAME?=lib/libmg_test.a
LIB_TEST_NAME:=$(LIB_TEST_NAME)

default: lib all

# Make dependecy list
dep:
	sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90
	$(MAKE) -C lib dep

lib:
	$(MAKE) -C lib lib

# Create list of all tests
TESTS=test1 test2 test3 test4 test5 test6 test7
.PHONY: $(TESTS)

all: lib $(TESTS)

# Ensure the library dependency
$(TESTS:=.o): | lib

test1: test1.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIB_TEST_NAME) $(LIB_NAME) $(LIBS)
ifdef RUN_TEST
	./$@
endif

test2: test2.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIB_TEST_NAME) $(LIB_NAME) $(LIBS)
ifdef RUN_TEST
	./$@
endif

test3: test3.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIB_TEST_NAME) $(LIB_NAME) $(LIBS)
ifdef RUN_TEST
	./$@
endif

test4: test4.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIB_TEST_NAME) $(LIB_NAME) $(LIBS)
ifdef RUN_TEST
	./$@
endif

test5: test5.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIB_TEST_NAME) $(LIB_NAME) $(LIBS)
ifdef RUN_TEST
	./$@
endif

test6: test6.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIB_TEST_NAME) $(LIB_NAME) $(LIBS)
ifdef RUN_TEST
	./$@
endif

test7: test7.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIB_TEST_NAME) $(LIB_NAME) $(LIBS)
ifdef RUN_TEST
	./$@
endif

clean:
	rm -f $(TESTS) $(TESTS:=.o) *.mod *.MOD *.cube *.VMG *.VMGASC *.nc
	$(MAKE) -C lib clean

# DO NOT DELETE THIS LINE - used by make depend
