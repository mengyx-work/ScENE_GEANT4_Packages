name := ScENE_June2013_10keVNeutron_Simplified



G4TARGET := $(name)
G4EXLIB := true

#G4WORKDIR :=.

#ifndef G4INSTALL
#  G4INSTALL = ../../..
#endif

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)
ROOTGLIBS       = $(shell root-config --glibs)

EXTRALIBS +=$(ROOTLIBS)
EXTRALIBS +=$(ROOTGLIBS)
#EXTRALIBS +=-L/usr/lib/CLHEP

CPPFLAGS += $(ROOTCFLAGS)

.PHONY: all
all:  lib bin

include $(G4INSTALL)/config/binmake.gmk

