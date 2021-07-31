bin_PROGRAMS = test

test_SOURCES = lion_class.F90  \
               intrinsic_class.F90  \
	       main.F90        

include_HEADERS = $(top_srcdir)/include/arrayofptrs_header.inc  \
		  $(top_srcdir)/include/cdlist_body.inc 	     \
		  $(top_srcdir)/include/cdlist_header.inc	     \
		  $(top_srcdir)/include/dlist_body.inc	     \
		  $(top_srcdir)/include/dlist_header.inc 	     \
		  $(top_srcdir)/include/list_body.inc	             \
		  $(top_srcdir)/include/list_header.inc	     \
		  $(top_srcdir)/include/prep.inc		     \
		  $(top_srcdir)/include/vector_body.inc	     \
		  $(top_srcdir)/include/vector_header.inc

AM_FCFLAGS = -I$(top_srcdir)/include/

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

all: test

# Propios
main.o: lion_class.o intrinsic_class.o


