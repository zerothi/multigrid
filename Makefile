.SUFFIXES: 
.SUFFIXES: .f90 .F90 .o .a
# We don't accept .f, nor .F files

# The phony list
.PHONY: default dep lib src test clean

default: src

# Make dependecy list
dep:
	@if [ -d src ] ;  then (cd src ; $(MAKE) dep) ; fi
	@if [ -d test ] ; then (cd test ; $(MAKE) dep) ; fi

lib: src
src: 
	@if [ -d src ] ;  then (cd src ; $(MAKE)) ; fi

test: lib
	@if [ -d test ] ;  then (cd test ; $(MAKE)) ; fi

clean:
	@if [ -d src ] ;  then (cd src ; $(MAKE) clean) ; fi
	@if [ -d test ] ; then (cd test ; $(MAKE) clean) ; fi
