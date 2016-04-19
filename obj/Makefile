PROJECT_ROOT = ..
include $(PROJECT_ROOT)/make_include
VPATH = $(SRC)

SOURCES := $(shell find $(SRC) -name '*.f90')
SOURCESPP := $(shell find $(SRC) -name '*.F90')
OBJS_ = $(SOURCES:%.f90=%.o) $(SOURCESPP:%.F90=%.o)
OBJS = $(OBJS_:$(SRC)/%=%)

all: init $(OBJS)

$(EXE):
	make -C $(PROJECT_ROOT) $(EXE)

ambient.o: ambient.f90
basal_surface.o: basal_surface.f90
cryosphere.o: cryosphere.f90 basal_surface.o glacier.o
glacier.o: glacier.f90
ground.o: ground.f90 basal_surface.o
ice_sheet.o: ice_sheet.f90 glacier.o
ice_shelf.o: ice_shelf.f90 glacier.o
plume.o: plume.f90 basal_surface.o

%.o: %.f90
	$(F90) $(FCFLAGS) -o $(ODIR)/$@ -c $<

%.o: %.F90
	$(F90) $(FCFLAGS) -o $(ODIR)/$@ -c $<

gcov:
	gcov *.gcno -pb

clean:
	rm -f *.gcov *.gcda *.gcno
	rm -f *~
	rm -f \#*
	rm -f *.o
	rm -f *.mod
