# LednickyEqn makefile

LIBS   = $(shell root-config --libs)
CFLAGS = $(shell root-config --cflags)
#IFLAGS = -I$(ROOTSYS)/include
#.PHONY: clean

all: main

main: lednickyEqn faddeeva lednickyInfo constraint pairsystem Main.cxx
	g++ Main.cxx Faddeeva.o LednickyEqn.o LednickyInfo.o ParameterConstraint.o PairSystem.o -o runMe $(LIBS) $(CFLAGS)

faddeeva: Faddeeva.cc
	g++ -c Faddeeva.cc

lednickyEqn: LednickyEqn.cxx 
	g++ -c LednickyEqn.cxx $(CFLAGS) 

lednickyInfo: LednickyInfo.cxx 
	g++ -c LednickyInfo.cxx $(CFLAGS) 

constraint: ParameterConstraint.cxx
	g++ -c ParameterConstraint.cxx $(CFLAGS)

pairsystem: PairSystem.cxx
	g++ -c PairSystem.cxx $(CFLAGS)



clean: 
	rm -f runMe LednickyEqn.o Faddeeva.o LednickyInfo.o PairSystem.o ParameterConstraint.o
