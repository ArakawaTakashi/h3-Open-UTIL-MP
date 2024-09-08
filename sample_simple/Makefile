include ../h3ou/sysdep/Makedef.$(SYSTEM)
INCCOUPLER = ../h3ou/include
LIBCOUPLER = ../h3ou/lib/libh3ou.a ../h3ou/lib/libjcup4.a
LIBCOUPLERC = ../h3ou/lib/libh3ouc.a ../h3ou/lib/libh3openp.a $(LIBCOUPLER)
all: compA compB compC

compA: compA.o
	$(LD) $(LDFLAGS) -o $@ compA.o $(LIBCOUPLER) 
compB: compB.o
	$(LD) $(LDFLAGS) -o $@ compB.o $(LIBCOUPLER) 
compC: compC.o
	$(LD) -nofor-main $(LDFLAGS) -o $@ compC.o $(LIBCOUPLERC) 


compA.o : compA.f90
	$(F90) $(CPPLOCAL) $(F90FLAGS) -I$(INCCOUPLER) -c compA.f90

compB.o : compB.f90
	$(F90) $(CPPLOCAL) $(F90FLAGS) -I$(INCCOUPLER) -c compB.f90

compC.o : compC.c
	$(CC) $(CPPLOCAL) $(CFLAGS) -I$(INCCOUPLER) -c compC.c

clean:
	rm *o ; rm *mod ; rm compA ; rm compB ; rm compC
