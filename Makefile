# $Id: Makefile 285 2012-01-02 16:02:07Z jbao $

SHELL = /bin/sh

SVNREV = $(shell svnversion -n)

CXX = $(HOME)/tool/openmpi-1.3.2/bin/mpic++
CC = gcc
CCC = g++
CFLAGS = -g 

#ifneq (, $(findstring geronto, $(HOST)))
#LIBS = -L/home/jbao/tool/openmpi-1.3.2/lib
#INCLUDES = -I/home/jbao/tool/openmpi-1.3.2/include
#endif

#LIBS = -lm
LIBS += -lgsl -lgslcblas -lm #-larmadillo -llapack -lblas 

OBJS = main.o dtw.o mds.o mscaling.o
TARGET = dtw mds hitmds2 interpolate pca corr

all: dtw info hitmds2 pca corr

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

dtw: main.o dtw.o
	$(CXX) $(CFLAGS) $(LIBS) $(INCLUDES) main.o dtw.o -o dtw

interpolate: interpolate.o
	$(CXX) $(CFLAGS) $(LIBS) $(INCLUDES) interpolate.o -o interpolate

mds: mds.o
	$(CXX) $(CFLAGS) $(LIBS) $(INCLUDES) mds.o -o mds

pca: pca.c
	$(CC) $(CFLAGS) $(INCLUDES) pca.c $(LIBS) -o pca

corr: mds.cpp
	$(CCC) $(CFLAGS) $(LIBS) $(INCLUDES) mds.cpp -o corr

info:
	@echo "Compiling. Stand by. 'not initialized' warnings are alright." && echo ""

# binary compilations
# 
hitmds2: hitmds2.c hitmds2.h
	$(CC) $(CFLAGS) -o hitmds2 hitmds2.c $(LIBS)

#$(TARGET): $(OBJS)
#	$(CXX) $(CFLAGS) $(LIBS) $(INCLUDES) $(OBJS) -o $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

dist:
	svn export . "dtw+mds.$(SVNREV)/"
	tar cvfz "dtw+mds.r$(SVNREV).tar.gz" "dtw+mds.$(SVNREV)/"
