# Structured Makefile for FindRoot
#
#

# Compiler and linker settings
FC = gfortran

# If you have local header files, list them here. (System headers are handled automatically.)
TARGET = Radius95
OBJS = func.o newt.o sec.o Radius95.o

# Debug flag: set debug=y on the command line to compile with debugging symbols.
debug ?= n
ifeq ($(debug), y)
    FFLAGS += -std=f2008 -O0 -g -Wall -Wextra -fimplicit-none -fcheck=all -fbacktrace
else
    FFLAGS += -std=f2008 -O2 -Wall -Wextra -fimplicit-none
endif

.PHONY: all clean pristine tags run

# Default target builds the executable and generates tags.
all: $(TARGET)

# Link the object file(s) into the executable "FindRoot"
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(TARGET)

#Compiles func.f90 into func.o
func.o: func.f90 
	$(FC) $(FFLAGS) -J. -I. -c func.f90

#Compiles newton.f90 into newton.o
newt.o: newt.f90 func.mod
	$(FC) $(FFLAGS) -J. -I. -c newt.f90

#Compiles sec.f90 into sec.o
sec.o: sec.f90 func.mod
	$(FC) $(FFLAGS) -J. -I. -c sec.f90

# Compiles main.f90 into main.o
Radius95.o: Radius95.f90 func.mod newt.mod sec.mod
	$(FC) $(FFLAGS) -J. -I. -c Radius95.f90

#Tell make that .mod files come from their .o files
func.mod: func.o; @true
newt.mod: newt.o; @true
sec.mod: sec.o; @true

#Runs the Program
run: $(TARGET)
	./$(TARGET)

# Remove object files.
clean:
	rm -f *.o *.mod

# Remove object files and the executable; then touch all files.
pristine: clean
	rm -f $(TARGET)
	touch *
	

# Generate tags for source and shader files.
tags:
	ctags *.f90



