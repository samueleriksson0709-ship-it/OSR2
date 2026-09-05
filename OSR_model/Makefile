# Cibles :
#   make           : compile tous les exécutables
#   make demo      : génère des données de test puis lance le pipeline
#   make test_osr  : compile et lance le test uniaxial uniquement
#   make clean     : nettoie les binaires et les fichiers .mod

FC = gfortran

# -ffree-line-length-none : several sources have lines past the 132 column limit
# of free form, which gfortran rejects by default (-Werror=line-truncation).
FFLAGS = -O3 -Wall -Wextra -fopenmp -ffree-line-length-none
LDFLAGS = -fopenmp

# Debug flags, used with: make debug
DEBUG_FLAGS = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace -fopenmp -ffree-line-length-none

MAIN_EXE = osr_pipeline_multi.exe

OBJS = osr_model.o pipeline_aux.o osr_pipeline_multi.o

# Default target

all: $(MAIN_EXE)

$(MAIN_EXE): $(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

# Object files

osr_model.o: osr_model.f90
		$(FC) $(FFLAGS) -c $<

pipeline_aux.o: pipeline_aux.f90 osr_model.o
		$(FC) $(FFLAGS) -c $<

osr_pipeline_multi.o: osr_pipeline_multi.f90 osr_model.o pipeline_aux.o
		$(FC) $(FFLAGS) -c $<

# Debug build

debug:
		$(MAKE) clean
		$(MAKE) FFLAGS="$(DEBUG_FLAGS)" LDFLAGS="-fopenmp"

# Utility targets

# rm on MSYS2 / Linux, del on a bare cmd.exe
clean:
		rm -f *.o *.mod *.exe || del /Q *.o *.mod *.exe

run: $(MAIN_EXE)
		.\$(MAIN_EXE)

.PHONY: all debug clean run
