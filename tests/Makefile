# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ..
include ../include.mk

# too lazy / dumb to go object by object in unit test deps
hal2sgOjbectsAll = $(wildcard ${rootPath}/*.o)
hal2sgOjbects=$(subst ../hal2sg.o,,${hal2sgOjbectsAll})
halTestSource = ${halRootPath}/api/tests/halAlignment*Test.cpp
halTestLib = ${halRootPath}/api/tests/halAlignment*Test.h
halTestInc = -I${halRootPath}/api/tests

all : unitTests

clean :
	rm -f *.o unitTests

unitTests : *.h *.cpp ${basicLibsDependencies} ../*.h ../*.cpp ${halTestLib} ${halTestSource}
	${cpp} -I../ ${cppflags} ${halTestInc} ${basicLibs} ${halTestSource} ${hal2sgOjbects} *.cpp -o unitTests

