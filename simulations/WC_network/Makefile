### Makefile for DDE_Example1 ###
PROGRAM_NAME := DDE_WC

# Dependencies
DEPEND := DDE_WC1_Driver.cpp DDE_WC1.cpp

# Set path to NFESOLVE directory (Change to wherever you have it installed)
NFESOLVE_DIR := $(HOME)/NFESOLVE-master

NFESOLVE_INCDIR := $(NFESOLVE_DIR)/include
NFESOLVE_LIBDIR := $(NFESOLVE_DIR)/lib

# Set up
CC := g++
CFLAGS := -std=c++11 -O2 -static-libgcc -static-libstdc++
OMP := -fopenmp
INC := -I $(NFESOLVE_INCDIR)
LIB_INCLUDES := -larmadillo


all: $(PROGRAM_NAME)

$(PROGRAM_NAME) : $(DEPEND)
	$(CC) $(CFLAGS) $(INC) -o $(PROGRAM_NAME) $^ -L$(NFESOLVE_LIBDIR) -lNFESOLVE $(LIB_INCLUDES)
	
par: $(PROGRAM_NAME)_PAR

$(PROGRAM_NAME)_PAR : $(DEPEND)
	$(CC) $(CFLAGS) $(OMP) $(INC) -o $(PROGRAM_NAME) $^ -L$(NFESOLVE_LIBDIR) -lNFESOLVE_PAR $(LIB_INCLUDES)


clean :
	rm -f $(PROGRAM_NAME)
