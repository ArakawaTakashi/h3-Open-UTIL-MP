H3OUHOME = /work/jh210022a/share/h3-Open-UTIL-MP/h3ou
include $(H3OUHOME)/sysdep/Makedef.$(SYSTEM)
INCCOUPLER = $(H3OUHOME)/include/$(SYSTEM)
LIBCOUPLER = $(H3OUHOME)/lib/$(SYSTEM)/libh3ou.a $(H3OUHOME)/lib/$(SYSTEM)/libjcup4.a
LIBCOUPLERC = $(H3OUHOME)/lib/$(SYSTEM)/libh3ouc.a $(H3OUHOME)/lib/$(SYSTEM)/libh3openp.a $(LIBCOUPLER)
all: compA compB compC

compA: compA.o
	$(LD) $(LDFLAGS) -o $@ compA.o $(LIBCOUPLER) 
	cp $@ $@.$(SYSTEM)
compB: compB.o
	$(LD) $(LDFLAGS) -o $@ compB.o $(LIBCOUPLER) 
	cp $@ $@.$(SYSTEM)
compC: compC.o
	$(LD) $(LDFLAGS) -o $@ compC.o $(LIBCOUPLERC) 
	cp $@ $@.$(SYSTEM)


compA.o : compA.f90
	$(F90) $(CPPLOCAL) $(F90FLAGS) -I$(INCCOUPLER) -c compA.f90

compB.o : compB.f90
	$(F90) $(CPPLOCAL) $(F90FLAGS) -I$(INCCOUPLER) -c compB.f90

compC.o : compC.c
	$(CC) $(CPPLOCAL) $(CFLAGS) -I$(INCCOUPLER) -c compC.c

clean:
	rm *o ; rm *mod ; rm compA ; rm compB ; rm compC
