include ./include.mk

sidegraphInc = sidegraph.h sgcommon.h sgsequence.h sgposition.h sgside.h sgjoin.h sgsegment.h

all : sgExport.a 

clean : 
	rm -f  sidegraph.o sglookup.o md5.o sgsql.o sgExport.a

unitTests : sgExport.a
	cd tests && make

sidegraph.o : sidegraph.cpp ${sidegraphInc}
	${cpp} ${cppflags} -I . sidegraph.cpp -c

sglookup.o : sglookup.cpp sglookup.h ${sidegraphInc}
	${cpp} ${cppflags} -I . sglookup.cpp -c

md5.o : md5.cpp md5.h
	${cpp} ${cppflags} -I . md5.cpp -c

sgsql.o : sgsql.cpp md5.h sgsql.h sglookup.h ${sidegraphInc}
	${cpp} ${cppflags} -I . sgsql.cpp -c

sgExport.a : sidegraph.o sglookup.o md5.o sgsql.o
	ar rc sgExport.a sidegraph.o sglookup.o md5.o sgsql.o

test : unitTests
	tests/unitTests
