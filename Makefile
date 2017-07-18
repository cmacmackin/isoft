#  
#  Copyright 2016 Christopher MacMackin <cmacmackin@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

# Directories
EXEC_FILE := main.o
SDIR := ./src
TDIR := ./tests
MDIR := ./mod
FDIR := ./factual
PFUNIT = $(HOME)/Code/pfunit-$(F90)

# Output
EXEC := isoft
TEXEC := tests.x
LIB := libisoft.a
FLIB := $(FDIR)/libfactual.a
FMDIR := $(FDIR)/mod
FMOD := $(FMDIR)/factual_mod.mod
PREFIX := $(HOME)/.local
LIBDIR := $(PREFIX)/lib
INCDIR := $(PREFIX)/include 
BINDIR := $(PREFIX)/bin

BLAS := openblas

# Include paths internal to project
PROJECT_INCDIRS := $(FMDIR) $(TDIR)

# The compiler
VENDOR ?= GNU
VENDOR_ = $(shell echo $(VENDOR) | tr A-Z a-z)

ifeq ($(VENDOR_),intel)
  F90      := ifort
  #FCFLAGS  := -O0 -g -DDEBUG -fpp -module $(MDIR) -traceback -assume realloc_lhs
  #LDFLAGS  := -O0 -g -traceback
  FCFLAGS  := -Ofast -fpp -module $(MDIR) -assume realloc_lhs
  LDFLAGS  := -Ofast
  COVFLAGS := 
else
  F90      := gfortran-6
  #FCFLAGS  := -O0 -g -pg -fcheck=all -DDEBUG -cpp -J$(MDIR) -ffpe-trap=overflow,invalid,zero
  #LDFLAGS  := -O0 -g -pg
  #COVFLAGS := -fprofile-arcs -ftest-coverage 
  FCFLAGS  := -O3 -cpp -J$(MDIR) -g
  LDFLAGS  := -O3
  COVFLAGS := 
endif

# Include paths
FCFLAGS += -I$(MDIR) $(PROJECT_INCDIRS:%=-I%) -I$(INCDIR) -I/usr/include -I$(PFUNIT)/mod

# Libraries for use at link-time
LIBS := -L$(LIBDIR) -lfftw3 -lnitsol -llapack95 -lhdf5_fortran -lhdf5hl_fortran -lflogging \
        -l$(BLAS) -llapack -lpenf
LDFLAGS += $(LIBS)

# A regular expression for names of modules provided by external libraries
# and which won't be contained in the module directory of this codebase
EXTERNAL_MODS := ^iso_(fortran_env|c_binding)|ieee_(exceptions|arithmetic|features)|openacc|omp_lib(_kinds)?|mpi|pfunit_mod|factual_mod|chebyshev_mod|array_pointer_mod|h5lt|hdf5|f95_lapack|logger_mod|penf|utils_mod$$

# Extensions of Fortran files, case insensitive
F_EXT := f for fpp f90 f95 f03 f08 f15

# Temporary work variables
_F_EXT := $(F_EXT) $(shell echo $(F_EXT) | tr a-z A-Z)
null :=
space := $(null) $(null)
EXT_PATTERN_GREP := '.*\.\($(subst $(space),\|,$(_F_EXT))\)'
EXT_PATTERN_SED := 's/([^ ]*)\.($(subst $(space),|,$(_F_EXT)))/\1.o/g;'

# Objects to compile
SRCS := $(shell find $(SDIR) -iregex $(EXT_PATTERN_GREP))
OBJS := $(shell find $(SDIR) -iregex $(EXT_PATTERN_GREP) | sed -r $(EXT_PATTERN_SED))
TOBJS := $(patsubst %.pf,%.o,$(wildcard $(TDIR)/*.pf))

META_FILE := $(SDIR)/meta_implementation.F90
META_OBJ  := $(META_FILE:$(suffix $(META_FILE))=.o)

#OBJS := $(filter-out $(META_OBJ), $(OBJS))

# "make" builds all
all: exec

all_objects: clean_src_mod $(OBJS)

exec: $(EXEC) 

$(EXEC): $(EXEC_FILE) $(LIB) $(FLIB)
	$(F90) $^ $(LDFLAGS) -o $@

install_exec: exec
	cp $(EXEC) $(BINDIR)

tests: $(TEXEC)
	./$(TEXEC)

$(TEXEC): $(TOBJS) $(LIB) $(FLIB) # $(TDIR)/testSuites.inc
	$(F90) -I$(PFUNIT)/include $(PFUNIT)/include/driver.F90 $(COVFLAGS) \
		$(FCFLAGS) $^ $(LIBS) -L$(PFUNIT)/lib -lpfunit -o $@

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) rcs $@ $(OBJS)

install_lib: lib
	cp $(LIB) $(LIBDIR)
	cp $(MDIR)/*.mod $(INCDIR)

$(FMOD): $(FLIB)

$(FLIB):
	make -C $(FLIB) all

ifeq ($(MAKECMDGOALS),clean)
else ifeq ($(MAKECMDGOALS),doc)
else
-include $(OBJS:.o=.d) $(TOBJS:.o=.d) $(EXEC_FILE:.o=.d)
endif

%.d: %.pf get_deps
	./get_deps $< $$\(MDIR\) "$(EXTERNAL_MODS)" $(PROJECT_INCDIRS) > $@

# General rule for building Fortran files, where $(1) is the file extension
define fortran_rule
%.o: %.$(1) $(FMOD) | $(MDIR)
	$$(F90) $$(COVFLAGS) $$(FCFLAGS) -c $$< -o $$@

%.d: %.$(1) get_deps
	./get_deps $$< \$$$$\(MDIR\) "$$(EXTERNAL_MODS)" $$(PROJECT_INCDIRS) > $$@
endef

# Register compilation rules for each Fortran extension
$(foreach EXT,$(_F_EXT),$(eval $(call fortran_rule,$(EXT))))

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $<  $@

$(META_OBJ): $(SRCS)

$(MDIR):
	@mkdir -p $@

clean: clean_obj clean_mod clean_deps clean_backups

clean_obj:
	@echo Deleting all object files
	@/bin/rm -rf $(OBJS) $(TOBJS)

clean_mod:
	@echo Deleting all module files
	@/bin/rm -rf $(MDIR)

clean_deps:
	@echo Deleting all dependency files
	@/bin/rm -rf $(OBJS:.o=.d) $(TOBJS:.o=.d)

clean_exec:
	@echo Deleting executable file
	@/bin/rm -rf $(EXEC)

clean_backups:
	@echo Deleting emacs backup files
	@/bin/rm -rf $(shell find '.' -name '*~') $(shell find '.' -name '\#*\#')

clean_src_mod:
	/bin/rm -rf $(SDIR)/*.mod

doc: documentation.md
	ford $<

.PHONEY: all all_objects exec install_exec tests lib install_lib clean \
	 clean_obj clean_mod clean_backups doc
