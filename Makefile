LOCAL_ROOT ?= $(HOME)/Code/local

VENDOR ?= GNU
VENDOR_ = $(shell echo $(VENDOR) | tr A-Z a-z)

ifeq ($(VENDOR_),gnu)
F90 = gfortran
PFUNIT = $(HOME)/Code/pfunit-gfortran
FCFLAGS = -O0 -g -J$(INC) -I$(SYSINC) -fcheck=all -DDEBUG
LDFLAGS = -O0 -g 
COVFLAGS = -fprofile-arcs -ftest-coverage
else ifeq ($(VENDOR_),intel)
F90 = ifort
PFUNIT = $(HOME)/Code/pfunit-ifort
FCFLAGS = -O0 -g -I$(INC) -module $(INC) -I$(SYSINC) -traceback \
          -assume realloc_lhs
LDFLAGS = -O0 -g -traceback
COVFLAGS = 
endif

#PFUNIT = /opt/pfunit-gfortran-5
export PFUNIT
export F90
export FCFLAGS
export LDFLAGS
export COVFLAGS

PWD = $(shell pwd)
INC = $(PWD)/mod
SYSINC = $(LOCAL_ROOT)/include
LIBS = -L$(PFUNIT)/lib -lpfunit $(COVFLAGS) -L. -lfactual -L$(LOCAL_ROOT)/lib -lfftw3 -lm


ARCHIVE = libfactual.a
SRC = ./src
TEST = ./tests
EXE = tests.x
#INC = ./mod
TESTINC = $(TEST)/mod

export ARCHIVE

# "make" builds all
all: lib

$(ARCHIVE): lib

lib:
	make -C $(SRC) lib
	cp $(SRC)/$(ARCHIVE) .

tests: $(EXE)
	./$(EXE) 

SUT: lib
	make -C $(TEST) tests

gcov: $(EXE)
	make -C $(SRC) gcov

%: $(ODIR)/%.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.F
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.f95
	$(FC) $(FCFLAGS) -o $@ -c $<

%.o: %.F95
	$(FC) $(FCFLAGS) -o $@ -c $<

clean:
	make -C $(SRC) clean
	make -C $(TEST) clean
	rm -f $(ARCHIVE) $(EXE)
	rm -f *.gcov *.gcda *.gcno

init:
	mkdir -p $(ODIR)
	mkdir -p $(INC)

echo:
	echo $(LOCAL_ROOT)
	echo $(PFUNIT)
	make -C $(SRC) echo


$(EXE): SUT
	$(F90) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/include \
		-I$(TEST) -I$(TESTINC) \
		$(PFUNIT)/include/driver.F90 $(TEST)/*.o $(FCFLAGS) \
		$(LIBS)

