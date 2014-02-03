.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

include ../arch.make

# The phony list
.PHONY: default deb clean

# We need to add the src-directory to the sources
FFLAGS += -I../src

LIB_NAME?=libmg.a
LIB_NAME:=../src/$(LIB_NAME)
default: all

# Make dependecy list
dep:
	sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90

# Grid objects
OBJS = m_cube.o

LIB_TEST_NAME?=libmg_test.a
$(LIB_TEST_NAME): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB_TEST_NAME) $(OBJS)
	-$(RANLIB) $(LIB_TEST_NAME)


all: $(LIB_TEST_NAME)
	for d in */ ; do \
		(cd $$d ; $(MAKE)) ; \
	done

clean:
	rm -f $(OBJS) *.mod *.MOD

# DO NOT DELETE THIS LINE - used by make depend
