#include make.inc


all: 
	make -C src
	@echo "done making the COECXS libraries"
	make -C examples
	@echo "done making the examples"

#install:
#	mkdir -p $(PREFIX)/include
#	cp include/* $(PREFIX)/include/.
#	mkdir -p $(PREFIX)/lib
#	cp lib/* $(PREFIX)/lib/.


clean: 
	make clean -C src
	make clean -C examples
	rm -f *~

clobber: 
	rm -f lib/libCOECXS.so lib/libCOECXS.a
	rm -f include/*.h


