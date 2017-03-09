.PHONY: clean

# set to the directory where Matlab is installed on your computer
MATLABROOT?=/usr/local/matlab

all:
	@echo
	@echo "	make linux32	to build on linux 32 bit"
	@echo "	make linux64	to build on linux 64 bit"
	@echo "	make macosx	to build on Mac OS X"
	@echo "	make clean	to erase unnecessary files"
	@echo "	make realclean	to erase all mex executables"
	@echo

linux32: clsolve.mexglx

linux64: clsolve.mexa64

macosx: clsolve.mexmaci64

tgz: drsolve.tgz

zip: drsolve.zip

clsolve.mexglx: clsolve.c fort.h fort.c
	mex clsolve.c fort.c -lmwblas -lmwlapack
	cp clsolve.mexglx ..

clsolve.mexa64: clsolve.c fort.h fort.c
	mex clsolve.c fort.c -lmwblas -lmwlapack
	cp clsolve.mexa64 ..

clsolve.mexmaci64: clsolve.c fort.h fort.c
	mex clsolve.c fort.c -lmwblas -lmwlapack
	cp clsolve.mexmaci64 ..

fort.h:
	cp ${MATLABROOT}/extern/examples/refbook/fort.h ./

fort.c:
	cp ${MATLABROOT}/extern/examples/refbook/fort.c ./

clean:
	rm -f *.mex*
	rm -f fort.h fort.c

realclean:
	rm -f *.mex*
	rm -f ../*.mex*
	rm -f fort.h fort.c

drsolve.tgz: realclean
	rm -f ../test/test*.mat
	(cd ../../; tar cvzf drsolve.tgz drsolve/)

drsolve.zip: realclean
	rm -f ../test/test*.mat
	(cd ../../; zip -r drsolve.zip drsolve/)

