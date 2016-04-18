PWD = $(shell pwd)
PROJECT_ROOT = .
include $(PROJECT_ROOT)/make_include

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

