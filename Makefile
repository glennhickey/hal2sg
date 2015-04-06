# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

all : hal2sg 

clean : 
	rm -f  hal2sg.o sidegraph.o sglookup.o hal2sg

hal2sg.o : hal2sg.cpp 
	${cpp} ${cppflags} -I . hal2sg.cpp -c

sidegraph.o : sidegraph.cpp sidegraph.h
	${cpp} ${cppflags} -I . sidegraph.cpp -c

sglookup.o : sidegraph.h sglookup.cpp sglookup.h
	${cpp} ${cppflags} -I . sglookup.cpp -c

hal2sg : hal2sg.o sidegraph.o sglookup.o ${basicLibsDependencies}
	${cpp} ${cppflags} ${basicLibs} hal2sg.o sidegraph.o sglookup.o -o hal2sg 

