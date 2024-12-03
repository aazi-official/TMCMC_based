SHELL = /bin/sh
TARGET = tmcmc_bimodal
CSRCS = main.cpp
OBJS = $(CSRCS:.c=.o)
DEPSRC = $(SRCS)
DEPCSRC = $(CSRCS)
LIBS = -L$(TMCMC_PATH)/lib -ltmcmc
INCL = -I$(TMCMC_PATH)/include
CPP = mpic++
CPPFLAGS = -O2 $(INCL)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CPP) $(CPPFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm -f *.o $(TARGET) *dat.* log_evidence.dat


