# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

all : hal2sg 

clean : 
	rm -f  hal2sg.o sidegraph.o sglookup.o hal2sg
	cd tests && make clean

unitTests : hal2sg
	cd tests && make

hal2sg.o : hal2sg.cpp 
	${cpp} ${cppflags} -I . hal2sg.cpp -c

sidegraph.o : sidegraph.cpp sidegraph.h sgcommon.h sgsequence.h sgposition.h sgside.h sgjoin.h
	${cpp} ${cppflags} -I . sidegraph.cpp -c

sglookup.o : sidegraph.o sglookup.cpp sglookup.h
	${cpp} ${cppflags} -I . sglookup.cpp -c

sgbuilder.o : sglookup.o sgbuilder.cpp sgbuilder.h ${basicLibsDependencies}
	${cpp} ${cppflags} -I . sgbuilder.cpp -c

hal2sg : hal2sg.o sidegraph.o sglookup.o sgbuilder.o ${basicLibsDependencies}
	${cpp} ${cppflags} ${basicLibs} hal2sg.o sidegraph.o sglookup.o sgbuilder.o -o hal2sg 

test : unitTests
	tests/unitTests

