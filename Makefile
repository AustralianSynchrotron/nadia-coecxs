#include make.inc


all: 
	make -C src
	@echo "done making the COECXS libraries"
	make -C tools
	@echo "done making the tools"
	make -C examples
	@echo "done making the examples"

#install:
#	mkdir -p $(PREFIX)/include
#	cp src/*.h $(PREFIX)/include
#	mkdir -p $(PREFIX)/lib
#	cp src/lib* $(PREFIX)/lib/
#	mkdir -p bin $(PREFIX)/bin/
#	cp tools/*.exe $(PREFIX)/bin/

clean: 
	make clean -C src
	make clean -C examples
	make clean -C tools
	rm -f *~

clobber:
	make clobber -C examples
	make clobber -C tools
	rm -f lib/libCOECXS.so lib/libCOECXS.a
	rm -f include/*.h


