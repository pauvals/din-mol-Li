# We are not following the standard layout of a GNU proyect
AUTOMAKE_OPTIONS = foreign


bin_PROGRAMS = dana

dana_SOURCES = dana.F90 \
               Constants.F90 \
               Cells.F90  \
               Neighbor.F90  \
               Elements.F90  \
               Program_Types.F90  \
               Errors.f90  \
               Strings.f90  \
               Groups.F90 

# Fortran Preprocesor Templates	       
FPTDDS=$(top_srcdir)/lib/fpt/include/
AM_FCFLAGS = -I$(FPTDDS)
include_HEADERS = $(FPTDDS)/prep.inc \
                  $(FPTDDS)/vector_header.inc \
                  $(FPTDDS)/vector_body.inc   \
                  $(FPTDDS)/list_header.inc   \
                  $(FPTDDS)/list_body.inc     \
                  $(FPTDDS)/dlist_header.inc  \
                  $(FPTDDS)/dlist_body.inc    \
                  $(FPTDDS)/cdlist_header.inc \
                  $(FPTDDS)/cdlist_body.inc   \
                  $(FPTDDS)/arrayofptrs_header.inc
 
DEBUG_FLAGS = -g -O0 -p
PROF_FLAGS = -p 
 
AM_FCFLAGS += $(OPENMP_CFLAGS) 
AM_FCFLAGS += -DDIM3

if COMPILER_IS_GFORTRAN
DEBUG_FLAGS += -fbacktrace
DEBUG_FLAGS += -fcheck=all
DEBUG_FLAGS += -Wall
DEBUG_FLAGS += -fbounds-check
DEBUG_FLAGS += -fcheck-array-temporaries 
DEBUG_FLAGS += -pedantic-errors
# F90DEBUG += -std=f2003 -pedantic 

# Some times needed for valgrind
# F90DEBUG += -static 

# To warn about use modules without only keyword
# F90DEBUG += -Wuse-without-only 
endif
    
if COMPILER_IS_IFORT
DEBUG_FLAGS += -check all 
DEBUG_FLAGS += -traceback 
DEBUG_FLAGS += -warn all  
DEBUG_FLAGS += -debug all 
#F90DEBUG += "-Mprof=func"
#F90DEBUG += "-Mprof=line"
endif

if DEBUG
  ENABLE_DEBUG=$(DEBUG_FLAGS)
endif

AM_FCFLAGS += $(ENABLE_DEBUG)



.PHONY: debug

debug:
	make ENABLE_DEBUG="$(DEBUG_FLAGS)"

all: dana

Constants.o:     
Strings.o:  Constants.o
Errors.o:   Constants.o Strings.o
Elements.o: Constants.o Strings.o Errors.o
 
Program_Types.o: Constants.o Errors.o Elements.o
dana.o: Groups.o Neighbor.o Errors.o
Program_Types.o: Constants.o Errors.o Elements.o
Groups.o:        Constants.o Elements.o Errors.o Program_Types.o
Neighbor.o:  Constants.o Groups.o Errors.o Cells.o
