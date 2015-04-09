# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ..
include ../include.mk

# too lazy / dumb to go object by object in unit test deps
hal2sgOjbectsAll = $(wildcard ${rootPath}/*.o)
hal2sgOjbects=$(subst ../hal2sg.o,,${hal2sgOjbectsAll})

all : unitTests

clean :
	rm -f *.o unitTests

unitTests : *.h *.cpp ${basicLibsDependencies} ../*.h ../*.cpp
	${cpp} -I../ ${cppflags} ${basicLibs} ${hal2sgOjbects} *.cpp -o unitTests

