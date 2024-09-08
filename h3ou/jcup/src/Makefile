TOPDIR = ../..
include ${TOPDIR}/Mkinclude

LIBJCUP = libjcup4.a

include files.mk

BLACK = c.black
WHITE = c.white

#MODS = $(JCUPOBJS:.o=.mod)
MODS = *.mod

all:
	cd impi ; $(MAKE) ; cd ..
	cd ${BLACK} ; $(MAKE) ; cd ..
	cd ${WHITE} ; $(MAKE) ; cd ..
	$(MAKE) $(JCUPLIBDIR)/$(LIBJCUP)
	$(MAKE) modules
	@echo "Complete making $(LIBJCUP)."

$(JCUPLIBDIR)/$(LIBJCUP): $(JCUPAPI)
#	$(AR) $(ARFLAGS) $@ $(ZOBJS) $(JCUPOBJS)
	$(AR) $(ARFLAGS) $@ $(JCUPOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

modules: $(JCUPAPI)
	$(INSTALL) $(MODS) $(JCUPINCDIR)


allclean: clean
	$(RM) -f $(LIBDIR)/$(LIBJCUP)

clean:
	cd impi ; make clean ; cd ..
	cd $(BLACK) ; make clean ; cd ..
	cd $(WHITE) ; make clean ; cd ..
	$(RM) -f *.o *.mod *.lst *.L



.SUFFIXES:
.SUFFIXES: .o .f90 .mod

.f90.o:
	$(FC) $(FFLAGS) -c $<

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)

