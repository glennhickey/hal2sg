# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ../
include ../include.mk

all : unitTests

clean :
	rm -f *.o unitTests

unitTests : ${rootPath}/hal2sg *.h *.cpp ${basicLibsDependencies}
	${cpp} ${cppflags} ${basicLibs} *.cpp -o unitTests

