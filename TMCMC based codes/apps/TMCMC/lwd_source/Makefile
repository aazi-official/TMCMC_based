SHELL = /bin/sh
TARGET = tmcmc_lwd
CSRCS = main.cpp
OBJS = $(CSRCS:.c=.o)
DEPSRC = $(SRCS)
DEPCSRC = $(CSRCS)

TMCMC_PATH = ../../../libs/TMCMC
FWLIB = ../../../libs/Forward
LIBS = -L$(TMCMC_PATH)/lib -ltmcmc_wb -L$(FWLIB)/lib -ldipole
INCL = -I$(TMCMC_PATH)/include -I$(FWLIB)/include

CPP = mpic++
CPPFLAGS = -O2 -ggdb -std=c++11 $(INCL)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CPP) $(CPPFLAGS) -o $(TARGET) $(OBJS) $(LIBS)
	mv $(TARGET) ../lwd_obj

clean:
	rm -f *.o $(TARGET) *dat.* log_evidence.dat



