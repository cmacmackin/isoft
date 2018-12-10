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

# System-specific variables
SYSTEM_LIB := /usr/lib      # Paths containing external libraries
SYSTEM_INC := /usr/include  # Paths containing module/include files
                            # for external libraries


SHELL := /bin/bash

# Directories
SDIR := src
TDIR := tests
MDIR := mod
LDIR := libs

# ISOFT path
IDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Included libraries
FACT_DIR := $(IDIR)$(LDIR)/factual
FACT_LIB := $(FACT_DIR)/libfactual.a
FACT_INC := $(FACT_DIR)/mod

FLOG_DIR := $(IDIR)$(LDIR)/flogging
FLOG_LIB := $(FLOG_DIR)/build/libflogging.a
FLOG_INC := $(FLOG_DIR)/build/include

FACE_DIR := $(IDIR)$(LDIR)/FACE
FACE_LIB := $(FACE_DIR)/lib/libface.a
FACE_INC := $(FACE_DIR)/lib/mod

LA_DIR := $(IDIR)$(LDIR)/lapack95
LA_LIB := $(LA_DIR)/lapack95.a
LA_INC := $(LA_DIR)/lapack95_modules 

NIT_DIR := $(IDIR)$(LDIR)/nitsol
NIT_LIB := $(NIT_DIR)/Nitsol/libnitsol.a
NIT_INC := 

PENF_DIR := $(IDIR)$(LDIR)/PENF
PENF_LIB := $(PENF_DIR)/static/penf.a
PENF_INC := $(PENF_DIR)/static/mod

PF_DIR := $(IDIR)$(LDIR)/pFUnit
PF_LIB := $(PF_DIR)/install/lib/libpfunit.a
PF_INC := $(PF_DIR)/install/mod
PF_DRIVE := $(PF_DIR)/install/include/driver.F90

# External libraries
BLAS := blas
LAPACK := lapack
HDF5 := hdf5_fortran
HDF5HL := hdf5hl_fortran
FFTW3 := fftw3

# Build tools
VENV := buildtools
FYPP := $(VENV)/bin/fypp
FOBIS := $(VENV)/bin/FoBiS.py
FORD := $(VENV)/bin/ford

# Output
TEXEC := tests.x
LIB := libisoft.a
SCRIPT_NAME := compile.sh

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
INC_FLAGS := -I$(FACT_INC) -I$(FLOG_INC) -I$(FACE_INC) \
	   -I$(LA_INC) -I$(PENF_INC) -I$(PF_INC) $(SYSTEM_INC:%=-I%)
FCFLAGS += -I$(MDIR) $(INC_FLAGS)

# Libraries for use at link-time
LIBS := $(FACT_LIB) $(FLOG_LIB) $(FACE_LIB) $(LA_LIB) $(NIT_LIB) $(PENF_LIB)
LIBDIRS := -L$(SYSTEM_LIB:%=-L%)
EXT_LIBS := -l$(FFTW3) -l$(HDF5) -l$(HDF5HL) -l$(BLAS) -l$(LAPACK)
LDFLAGS += $(LIBS) $(LIBDIRS) $(EXT_LIBS)

# Compile script
define COMP_SCRIPT
#!/bin/bash

if [ "$$#" -gt 0 ]
then
    ifile=$$1
else
    ifile="main.f90"
fi
if [ "$$#" -gt 1 ]
then
    iexec=$$2
else
    iexec="isoft"
fi

f90="$(F90)"
idir="$(IDIR)"
ilib="$${idir}/$(LIB)"

flags="-O3 -I$(IDIR)$(MDIR) $(INC_FLAGS) $(LDFLAGS)"

command="$${f90} $${ifile} -o $${iexec} $${ilib} $${flags}"
echo $$command
$$($$command)
endef
export COMP_SCRIPT

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
TOBJS := $(patsubst %.pf,%.o,$(wildcard $(TDIR)/*.pf)) $(patsubst %.f90,%.o,$(wildcard $(TDIR)/*.f90))

META_FILE := $(SDIR)/meta_implementation.F90
META_OBJ  := $(META_FILE:$(suffix $(META_FILE))=.o)

#OBJS := $(filter-out $(META_OBJ), $(OBJS))

# "make" builds all
all: lib

all_objects: clean_src_mod $(OBJS)

tests: $(TEXEC)
	./$(TEXEC)

$(TEXEC): $(TOBJS) $(LIB) $(PF_LIB) | $(TDIR)/testSuites.inc
	$(F90) -I$(TDIR) $(PF_DRIVE) $(COVFLAGS) \
		$(FCFLAGS) $^ $(LDFLAGS) -o $@

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) rcs $@ $(OBJS)

$(FACT_LIB): $(FYPP)
	source $(VENV)/bin/activate
	make -C $(FACT_DIR)

$(FACE_LIB): $(FOBIS)
	source $(VENV)/bin/activate
	cd $(FACE_DIR); \
	FoBiS.py build -f fobos -mode face-static-gnu

$(FLOG_LIB): $(FACE_LIB) $(FOBIS)
	source $(VENV)/bin/activate
	cd $(FLOG_DIR); \
	FoBiS.py build -f fobos -mode gnu-static -i ../../$(FACE_INC)

$(LA_LIB): | $(LA_INC)
	make -C $(LA_DIR)/SRC single_double_complex_dcomplex OPTS0='-O3'

$(LA_INC):
	mkdir $(LA_INC)

$(NIT_LIB):
	make -C $(NIT_DIR)

$(PENF_LIB): $(FOBIS)
	source $(VENV)/bin/activate
	cd $(PENF_DIR); \
	FoBiS.py build -f fobos -mode static-gnu

$(PF_LIB):
	make -C $(PF_DIR) tests F90=$(F90) F90_VENDOR=GNU
	make -C $(PF_DIR) install INSTALL_DIR=install

$(FOBIS): $(VENV)
	source $(VENV)/bin/activate
	pip install FoBiS.py

$(FORD): $(VENV)
	source $(VENV)/bin/activate
	pip install ford

$(FYPP): $(VENV)
	source $(VENV)/bin/activate
	pip install fypp

$(VENV):
	virtualenv $(VENV)

ifeq ($(MAKECMDGOALS),clean)
else ifeq ($(MAKECMDGOALS),doc)
else
-include $(OBJS:.o=.d) $(TOBJS:.o=.d) $(EXEC_FILE:.o=.d)
endif

%.d: %.pf get_deps
	./get_deps $< $$\(MDIR\) "$(EXTERNAL_MODS)" $(PROJECT_INCDIRS) > $@

# General rule for building Fortran files, where $(1) is the file extension
define fortran_rule
%.o: %.$(1) $(LIBS) | $(MDIR)
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

script: $(SCRIPT_NAME)

$(SCRIPT_NAME): | Makefile
	@echo "$$COMP_SCRIPT" > compile.sh
	@chmod a+x $@
	@echo Created script compile.sh

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

doc: documentation.md $(FORD)
	source $(VENV)/bin/activate
	ford $<

.PHONEY: all all_objects exec install_exec tests lib install_lib clean \
	 clean_obj clean_mod clean_backups doc
