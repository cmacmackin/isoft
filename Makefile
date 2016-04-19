PROJECT_ROOT = .
include $(PROJECT_ROOT)/make_include

all: $(EXE)

$(EXE): objs
	$(F90) -o $@ $(ODIR)/*.o $(LDFLAGS) $(LIBS) $(COVFLAGS)

libfactual:
	make -C $(FACTUAL) lib

objs: libfactual
	make -C $(ODIR) all

tests:  $(TESTEXE)
	$(TESTEXE) 

$(TESTEXE): SUT
	$(F90) -o $@ -I$(TEST) $(PFUNIT)/include/driver.F90 $(TESTODIR)/*.o \
		$(ODIR)/*.o $(LDFLAGS) $(LIBS) -L$(PFUNIT)/lib -lpfunit \
		$(COVFLAGS)

SUT:    $(EXE)
	make -C $(TEST) tests

gcov:   $(TESTEXE)
	make -C $(ODIR) gcov

clean:
	make -C $(FACTUAL) clean
	make -C $(ODIR) clean
	make -C $(SRC) clean
	make -C $(TEST) clean
	rm -f $(TESTEXE) $(EXE)
	rm -f *.gcov *.gcda *.gcno
	rm -f *~
	rm -f \#*
	rm -f $(INC)/*.mod
	rm -f $(ODIR)/*.o
