# $Id: Makefile,v 1.9 2014/10/17 12:15:01 heine Exp $
# Last modified on February 15, 2017 by M. Paun

# for mac terminal
#CXX=g++-6
#CXXFLAGS=-O2 -g -Wall -D_REENTRANT

#FC=gfortran-6
#FFLAGS=-O2 -g -Wall
#FLIBS=-lgfortran -lquadmath

######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
## IBAMR_SRC_DIR = /xlwork6/hgao/IBAMR2017/IBAMR-poroelasticity/IBAMR
## IBAMR_BUILD_DIR = /xlwork6/hgao/IBAMR2017/IBAMR-poroelasticity/opt

IBAMR_SRC_DIR = /xlwork6/hgao/IBAMR2019Feb/IBAMR-git/IBAMR
IBAMR_BUILD_DIR = /xlwork6/hgao/IBAMR2019Feb/IBAMR-git/opt

## IBAMR_SRC_DIR = /xlwork1/scott/sfw/ibamr/IBAMR
## IBAMR_BUILD_DIR = /xlwork1/scott/sfw/ibamr/ibamr-objs-opt


######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

#for windows terminal and different compiler
CC=mpicc
CXX=mpicxx
CXXFLAGS +=-O2 -g -Wall -D_REENTRANT -fPIC

ifeq ($(shell arch),ppc)
	FC=/usr/local/bin/gfortran
	FLIBS=-L$(shell $(FC) --print-file-name=) -lf95
	FFLAGS=-O2 -g -Wall -fno-second-underscore
else
	FC=/usr/local/bin/gfortran
	FLIBS=-lgfortran
	FFLAGS=-O2 -g -Wall -fPIC
endif

FC=gfortran



SRC = $(wildcard *.C) $(wildcard /home/staff1/nthekkethi/FEMLDE_7/*.C)
PDIM = 3
OBJS = $(SRC:%.C=%.o) $(IBAMR_LIB_3D) $(IBTK_LIB_3D)

main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) -o main

# OBJS1=impedance_sub.o new_match.o f90_tools.o external_sub.o ext_new_match.o

# main: $(OBJS) $(OBJS1)
# 	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(OBJS1) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) -o main
#
#
# new_match.o: new_match.f90 f90_tools.o
# 	$(FC) -c $(FFLAGS) new_match.f90
#
# f90_tools.o: f90_tools.f90
# 	$(FC) -c $(FFLAGS) f90_tools.f90
#
# impedance_sub.o: impedance_sub.f90 f90_tools.o new_match.o
# 	$(FC) -c $(FFLAGS) impedance_sub.f90
#
# ext_new_match.o: ext_new_match.f90 f90_tools.o
# 	$(FC) -c $(FFLAGS) ext_new_match.f90
#
# external_sub.o: external_sub.f90 f90_tools.o ext_new_match.o
# 	$(FC) -c $(FFLAGS) external_sub.f90

clean:
	$(RM) main
	$(RM) *.o *.lo *.objs *.ii *.int.c
	$(RM) -r .libs

-include $(SRC:%.C=%.d)
