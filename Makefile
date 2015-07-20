# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

sidegraphInc = ${sgExportPath}/sidegraph.h ${sgExportPath}/sgcommon.h ${sgExportPath}/sgsequence.h ${sgExportPath}/sgposition.h ${sgExportPath}/sgside.h ${sgExportPath}/sgjoin.h ${sgExportPath}/sgsegment.h

all : hal2sg 

clean : 
	rm -f  hal2sg.o sglookback.o snphandler.o sgbuilder.o halsgsql.o hal2sg
	cd sgExport && make clean
	cd tests && make clean

unitTests : hal2sg
	cd tests && make

hal2sg.o : hal2sg.cpp halsgsql.h sgbuilder.h ${sgExportPath}/sglookup.h snphandler.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . hal2sg.cpp -c

sglookback.o : sglookback.cpp sglookback.h ${sgExportPath}/sglookup.h ${sidegraphInc}
	${cpp} ${cppflags} -I . sglookback.cpp -c

snphandler.o : snphandler.cpp snphandler.h sglookback.h ${sgExportPath}/sglookup.h ${sidegraphInc}
	${cpp} ${cppflags} -I . snphandler.cpp -c

sgbuilder.o : sgbuilder.cpp sgbuilder.h ${sgExportPath}/sglookup.h sglookback.h snphandler.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . sgbuilder.cpp -c

halsgsql.o : halsgsql.cpp halsgsql.h ${sgExportPath}/sglookup.h ${sgExportPath}/sgsql.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . halsgsql.cpp -c

${sgExportPath}/sgExport.a : ${sgExportPath}/*.cpp ${sgExportPath}/*.h
	cd ${sgExportPath} && make

hal2sg :  hal2sg.o sglookback.o snphandler.o sgbuilder.o halsgsql.o ${basicLibsDependencies}
	${cpp} ${cppflags} ${basicLibs} hal2sg.o sglookback.o snphandler.o sgbuilder.o halsgsql.o -o hal2sg 

test : unitTests
	pushd .  && cd ${sgExportPath} && make test && popd && tests/unitTests

