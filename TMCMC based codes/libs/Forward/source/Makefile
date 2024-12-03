CC = g++
CXX = g++
CCFLAGS +=  -g -c -lm -std=c++0x -O0 -fPIC
OBJS:= forward
LIB = libdipole.a
.SUFFIXES : .o .cpp
INCDIR = ../include/
INCPATH = -I$(INCDIR)
all : $(LIB)
$(OBJS): Curves.o FastTIDriver.o QuadRule.o TxLine.o
	$(CC) Curves.o \
        FastTIDriver.o \
        QuadRule.o \
        TxLine.o \
         -fopenmp -c -o forward


TxLine.o : TxLine.cpp
	$(CXX) TxLine.cpp $(INCPATH) $(CCFLAGS)
QuadRule.o : QuadRule.cpp TxLine.o
	$(CXX) QuadRule.cpp TxLine.o $(INCPATH) $(CCFLAGS)
FastTIDriver.o : FastTIDriver.cpp TxLine.o QuadRule.o
	$(CXX) FastTIDriver.cpp TxLine.o QuadRule.o $(INCPATH) $(CCFLAGS)
Curves.o : Curves.cpp FastTIDriver.o 
	$(CXX) Curves.cpp FastTIDriver.o $(INCPATH) $(CCFLAGS)

$(LIB) : Curves.o FastTIDriver.o QuadRule.o TxLine.o
	ar cr $@ Curves.o FastTIDriver.o QuadRule.o TxLine.o
	mv $(LIB) ../lib/


clean :
	rm -f *~
	rm -f *.o
