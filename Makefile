# LednickyEqn makefile

IFLAGS = -I$(ROOTSYS)/include
LFLAGS = -L$(ROOTSYS)/lib -lHist -lMathCore -lCore -lRIO -lGpad -lMinuit
#.PHONY: clean

all: lednicky

lednicky: LednickyEqn.cxx 
	g++ LednickyEqn.cxx -o lednicky $(IFLAGS) $(LFLAGS)

