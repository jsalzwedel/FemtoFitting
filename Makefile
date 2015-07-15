# LednickyEqn makefile

LIBS   = $(shell root-config --libs)
CFLAGS = $(shell root-config --cflags)
#IFLAGS = -I$(ROOTSYS)/include
#.PHONY: clean

all: main

main: lednicky faddeeva
	g++ Faddeeva.o LednickyEqn.o -o lednicky $(LIBS)

faddeeva: Faddeeva.cc
	g++ -c Faddeeva.cc

lednicky: LednickyEqn.cxx 
	g++ -c LednickyEqn.cxx $(CFLAGS) 

clean: 
	rm -f lednicky Lednicky.o Faddeeva.o
