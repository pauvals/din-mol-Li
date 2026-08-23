bin_PROGRAMS = dana

dana_SOURCES = dana.F90 \
               gems_core/Constants.F90 \
               gems_core/Strings.f90  \
               gems_core/Errors.f90  \
               gems_core/Elements.F90  \
               gems_core/Math.F90  \
               gems_core/Program_Types.F90  \
               gems_core/Groups.F90 \
               gems_core/Halo.F90 \
               gems_core/Cells.F90  \
               gems_core/Neighbor.F90
    
# zeros...
dana_LDADD =
AM_FCFLAGS =

AM_FCFLAGS += $(OPENMP_CFLAGS)
AM_FCFLAGS += $(COARRAY_FCFLAGS)
AM_FCFLAGS += -DDIM3
dana_LDADD += $(COARRAY_LIBS)
        
# External software
# =================

# Note on local vs external libraries:
# > We recommend that you avoid using -l options in LDADD or prog_LDADD when
#   referring to libraries built by your package. Instead, write the file name of
#   the library explicitly as in the above cpio example. Use -l only to list
#   third-party libraries. If you follow this rule, the default value of
#   prog_DEPENDENCIES will list all your local libraries and omit the other ones.
# [source](https://www.gnu.org/software/automake/manual/html_node/Linking.html)

# Fortran Preprocesor Templates
# -----------------------------
AM_FCFLAGS += -I$(top_srcdir)/lib/fpt/include/

# Debug
# =====
   
DEBUG_FLAGS = -g -O0 -p
PROF_FLAGS = -p 

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
gems_core/Strings.o:  gems_core/Constants.o
gems_core/Errors.o:   gems_core/Constants.o gems_core/Strings.o
gems_core/Elements.o: gems_core/Constants.o gems_core/Strings.o gems_core/Errors.o
 
gems_core/Math.o: gems_core/Constants.o gems_core/Errors.o
gems_core/Program_Types.o: gems_core/Constants.o gems_core/Errors.o gems_core/Elements.o
gems_core/Groups.o:        gems_core/Constants.o gems_core/Elements.o gems_core/Errors.o gems_core/Math.o gems_core/Program_Types.o
gems_core/Halo.o:          gems_core/Constants.o gems_core/Groups.o gems_core/Program_Types.o
gems_core/Cells.o:         gems_core/Program_Types.o gems_core/Constants.o gems_core/Groups.o gems_core/Errors.o gems_core/Halo.o
gems_core/Neighbor.o:      gems_core/Program_Types.o gems_core/Constants.o gems_core/Groups.o gems_core/Errors.o gems_core/Cells.o gems_core/Halo.o
dana.o: gems_core/Groups.o gems_core/Neighbor.o gems_core/Errors.o gems_core/Elements.o gems_core/Constants.o

clean-local:
	rm -f *.mod
