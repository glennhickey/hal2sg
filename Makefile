# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

sidegraphInc = sidegraph.h sgcommon.h sgsequence.h sgposition.h sgside.h sgjoin.h sgsegment.h

all : hal2sg 

clean : 
	rm -f  hal2sg.o sidegraph.o sglookup.o snphandler.o sgbuilder.o md5.o sgsql.o hal2sg
	cd tests && make clean

unitTests : hal2sg
	cd tests && make

hal2sg.o : hal2sg.cpp sgsql.h sgbuilder.h sglookup.h snphandler.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . hal2sg.cpp -c

sidegraph.o : sidegraph.cpp ${sidegraphInc}
	${cpp} ${cppflags} -I . sidegraph.cpp -c

sglookup.o : sglookup.cpp sglookup.h ${sidegraphInc}
	${cpp} ${cppflags} -I . sglookup.cpp -c

snphandler.o : snphandler.cpp snphandler.h ${sidegraphInc}
	${cpp} ${cppflags} -I . snphandler.cpp -c

sgbuilder.o : sgbuilder.cpp sgbuilder.h sglookup.h snphandler.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . sgbuilder.cpp -c

md5.o : md5.cpp md5.h
	${cpp} ${cppflags} -I . md5.cpp -c

sgsql.o : sgsql.cpp md5.h sgsql.h sglookup.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . sgsql.cpp -c

hal2sg : hal2sg.o sidegraph.o sglookup.o snphandler.o sgbuilder.o md5.o sgsql.o ${basicLibsDependencies}
	${cpp} ${cppflags} ${basicLibs} hal2sg.o sidegraph.o sglookup.o snphandler.o sgbuilder.o md5.o sgsql.o -o hal2sg 

test : unitTests
	tests/unitTests

