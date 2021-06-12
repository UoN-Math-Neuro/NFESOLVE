CC := g++
CFLAGS := -std=c++11 -O2 -Wl,-subsystem,console
OMP := -fopenmp
INC := -I include/NFESOLVE/ $(addprefix -I , $(dir $(wildcard include/*/*/))) 
LIBS := -larmadillo 

SRCDIR := src
BUILDDIR_SER := build/ser
$(shell   mkdir -p $(BUILDDIR_SER))
BUILDDIR_PAR := build/par
$(shell   mkdir -p $(BUILDDIR_PAR))
LIBDIR := lib
$(shell   mkdir -p $(LIBDIR))

SRCEXT := cpp

SOURCES := $(wildcard $(SRCDIR)/**/*$(SRCEXT))
OBJECTS_SER := $(addprefix $(BUILDDIR_SER)/, $(notdir $(SOURCES:.$(SRCEXT)=.o)))
OBJECTS_PAR := $(addprefix $(BUILDDIR_PAR)/, $(notdir $(SOURCES:.$(SRCEXT)=.o)))

VPATH := $(wildcard $(SRCDIR)/*)

all : NFESOLVE NFESOLVE_PAR
.PHONY: all

$(BUILDDIR_SER)/%.o : %.$(SRCEXT)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@ $(LIBS)

NFESOLVE : $(OBJECTS_SER)
	ar rcs $(LIBDIR)/libNFESOLVE.a $(OBJECTS_SER)
	ranlib $(LIBDIR)/libNFESOLVE.a

$(BUILDDIR_PAR)/%.o : %.$(SRCEXT)
	$(CC) $(CFLAGS) $(OMP) $(INC) -c $< -o $@ $(LIBS)

NFESOLVE_PAR : $(OBJECTS_PAR)
	ar rcs $(LIBDIR)/libNFESOLVE_PAR.a $(OBJECTS_PAR)
	ranlib $(LIBDIR)/libNFESOLVE_PAR.a


cleanall : clean
	rm -f $(LIBDIR)/*

clean :
	@echo " Cleaning..."
	rm -f $(BUILDDIR_SER)/*
	rm -f $(BUILDDIR_PAR)/*

clean_ser :
	@echo " Cleaning..."
	rm -f $(BUILDDIR_SER)/*

clean_par :
	@echo " Cleaning..."
	rm -f $(BUILDDIR_PAR)/*

.PHONY: clean

remake : clean all
remake_ser: clean_ser NFESOLVE
remake_par: clean_par NFESOLVE_PAR
