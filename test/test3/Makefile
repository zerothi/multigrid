.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

include ../../arch.make

# The phony list
.PHONY: default deb clean

# Add include flags
FFLAGS += -I../
FFLAGS += -I../../src

# Define the libraries
LIB_NAME?=libmg.a
LIB_NAME:=../../src/$(LIB_NAME)
LIB_TEST_NAME?=libmg_test.a
LIB_TEST_NAME:=../$(LIB_TEST_NAME)

default: test1

# Make dependecy list
dep:
	sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90

test1: four.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIB_TEST_NAME) $(LIB_NAME)
	./$@

clean:
	rm -f test1 $(OBJS) *.mod *.MOD

# DO NOT DELETE THIS LINE - used by make depend
