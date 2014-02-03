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

default: test

# Make dependecy list
dep:
	sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90

# Grid objects
OBJS = topbottom.o

test: $(OBJS)
	$(FC) -o test $(OBJS) $(LIB_TEST_NAME) $(LIB_NAME)
	./test

clean:
	rm -f test $(OBJS) *.mod *.MOD

# DO NOT DELETE THIS LINE - used by make depend
