SHELL = /bin/sh
TARGET = rjTMCMC
CSRCS = rjTMCMC.cpp
OBJS = $(CSRCS:.c=.o)
DEPSRC = $(SRCS)
DEPCSRC = $(CSRCS)

TRJMCMCLIB = ../../../libs/T-RJMCMC
FWLIB = ../../../libs/Forward
LIBS = -L$(TRJMCMCLIB)/lib -ltrjmcmc -L$(FWLIB)/lib -ldipole
INCL = -I$(TRJMCMCLIB)/include -I$(FWLIB)/include

CPP = mpic++
CPPFLAGS = -O0 -ggdb -std=c++11 $(INCL)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CPP) $(CPPFLAGS) -o $(TARGET) $(OBJS) $(LIBS)
	mv $(TARGET) ../lwd_obj

clean:
	rm -f *.o ../lwd_obj/$(TARGET) *dat.* log_evidence.dat



