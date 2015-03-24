.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

include arch.make

# The phony list
.PHONY: default dep lib src ncdf test clean mg

default: lib

ifneq (,$(findstring -DCDF,$(FPPFLAGS)))
ncdf:
	$(MAKE) -C ncdf "ARCH_MAKE=$$(pwd)/arch.make"
else
ncdf:
	@echo "Compiling without NetCDF support"
endif



# Make dependency list
dep:
	$(MAKE) -C src dep
	$(MAKE) -C test dep

lib: src
src: ncdf
	$(MAKE) -C src

mg: ncdf
	$(MAKE) -C src mg

test: lib
	$(MAKE) -C test

clean:
	$(MAKE) -C src clean
	$(MAKE) -C test clean
	-$(MAKE) -C ncdf clean "ARCH_MAKE=$$(pwd)/arch.make" clean
