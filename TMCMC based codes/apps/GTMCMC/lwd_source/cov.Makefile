SHELL = /bin/sh
TARGET = gtmcmc_cov_lwd
CSRCS = main_cov.cpp
OBJS = $(CSRCS:.c=.o)
DEPSRC = $(SRCS)
DEPCSRC = $(CSRCS)

FWLIB = ../../../libs/Forward
GTMCMC_PATH = ../../../libs/GTMCMC
LIBS = -L$(GTMCMC_PATH)/lib -lgtmcmc_cov -L$(FWLIB)/lib -ldipole
INCL = -I$(GTMCMC_PATH)/include -I$(FWLIB)/include

CPP = mpic++
CPPFLAGS = -O2 -ggdb -std=c++11 $(INCL)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CPP) $(CPPFLAGS) -o $(TARGET) $(OBJS) $(LIBS)
	mv $(TARGET) ../lwd_obj

clean:
	rm -f *.o $(TARGET) *dat.* log_evidence.dat



