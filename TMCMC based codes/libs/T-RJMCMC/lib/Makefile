SHELL = /bin/sh
TARGET = libtrjmcmc.a
CSRCS = tmcmc.cpp
DEPCSRC = $(CSRCS)
OBJS = $(CSRCS:.cpp=.o)
INCFILES = tmcmc.h
INCDIR=../include
INCL = -I. 
CPP = mpic++
CPPFLAGS = -O0 -std=gnu++11 -ggdb $(INCL)
RANLIB = ranlib
AR = ar cr

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)


clean:
	rm -f $(OBJS) $(TARGET)

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $*.cpp
