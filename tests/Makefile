### Makefile for tests ###
PROGRAM_NAME := DelayUtilTests

# Dependencies
DEPEND := t_DelayUtils.cpp

# Set path to NFESOLVE directory (Change to wherever you have it installed)
NFESOLVE_DIR := $(HOME)/NFESOLVE-master/

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

clean :
	rm -f $(PROGRAM_NAME)
