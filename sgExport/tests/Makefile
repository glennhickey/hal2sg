rootPath = ..
include ${rootPath}/include.mk

all : unitTests

clean :
	rm -f *.o unitTests

unitTests : CuTest.o unitTests.o sidegraphTests.o sglookupTests.o 

CuTest.o : CuTest.h CuTest.c
	${cxx} ${cflags} -c CuTest.c	

unitTests.o : unitTests.h unitTests.cpp
	${cpp} ${cppflags} -I${rootPath}/ unitTests.cpp -c

sidegraphTests.o : sidegraphTests.cpp ${rootPath}/*.h
	${cpp} ${cppflags} -I${rootPath}/ sidegraphTests.cpp -c

sglookupTests.o : sglookupTests.cpp ${rootPath}/*.h
	${cpp} ${cppflags} -I${rootPath}/ sglookupTests.cpp -c

unitTests : CuTest.o unitTests.o sidegraphTests.o sglookupTests.o ${rootPath}/sgExport.a
	${cpp} -I${rootPath}/ ${cppflags} ${rootPath}/sgExport.a CuTest.o unitTests.o sidegraphTests.o sglookupTests.o -o unitTests

