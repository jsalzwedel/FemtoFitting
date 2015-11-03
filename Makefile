# LednickyEqn makefile

LIBS   = $(shell root-config --libs) -lMinuit
CFLAGS = $(shell root-config --cflags) -Wall -O2
#IFLAGS = -I$(ROOTSYS)/include
#.PHONY: clean

all: runMe

runMe: LednickyEqn.o Faddeeva.o LednickyInfo.o ParameterConstraint.o PairSystem.o Fitter.o Main.cxx Faddeeva.cc LednickyEqn.cxx LednickyInfo.cxx ParameterConstraint.cxx PairSystem.cxx Fitter.cxx
	g++ Main.cxx Faddeeva.o LednickyEqn.o LednickyInfo.o ParameterConstraint.o PairSystem.o Fitter.o -o runMe $(LIBS) $(CFLAGS) 

Faddeeva.o: Faddeeva.cc
	g++ -c Faddeeva.cc

LednickyEqn.o: LednickyEqn.cxx 
	g++ -c LednickyEqn.cxx $(CFLAGS)

LednickyInfo.o: LednickyInfo.cxx 
	g++ -c LednickyInfo.cxx $(CFLAGS) 

ParameterConstraint.o: ParameterConstraint.cxx
	g++ -c ParameterConstraint.cxx $(CFLAGS)

PairSystem.o: PairSystem.cxx
	g++ -c PairSystem.cxx $(CFLAGS)

Fitter.o: Fitter.cxx
	g++ -c Fitter.cxx $(CFLAGS)



clean: 
	rm -f runMe *.o
