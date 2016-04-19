PROJECT_ROOT = ..
include $(PROJECT_ROOT)/make_include
clean:
	rm -f *.gcov *.gcda *.gcno
	rm -f *~
	rm -f \#*
	rm -f *.o
	rm -f *.mod

%.o: %.f90
	make -C $(ODIR) $@

%.o: %.F90
	make -C $(ODIR) $@
