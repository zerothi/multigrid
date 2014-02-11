.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

include ../arch.make

# The phony list
TESTS=test1 test2 test3
TESTS=test1
.PHONY: default lib deb clean all $(TESTS)

# We need to add the src-directory to the sources
FFLAGS += -I../src

default: all

# Make dependecy list
dep:
	sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90
	for d in */ ; do \
		(cd $$d ; $(MAKE) dep) ; \
	done

# Grid objects
OBJS = m_cube.o m_mg_info.o m_time.o

LIB_TEST_NAME?=libmg_test.a
lib: $(LIB_TEST_NAME)
$(LIB_TEST_NAME): $(OBJS)
	$(AR) $(ARFLAGS) $(LIB_TEST_NAME) $(OBJS)
	-$(RANLIB) $(LIB_TEST_NAME)


all: $(LIB_TEST_NAME) $(TESTS)

test1:
	(cd $@ ; $(MAKE))

test2:
	(cd $@ ; $(MAKE))

test3:
	(cd $@ ; $(MAKE))

clean:
	rm -f $(OBJS) *.mod *.MOD
	for d in */ ; do \
		(cd $$d ; $(MAKE) clean) ; \
	done

# DO NOT DELETE THIS LINE - used by make depend
