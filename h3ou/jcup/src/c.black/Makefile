TOPDIR = ../../..
include ${TOPDIR}/Mkinclude
include files.mk

CMODS = $(COBJS:.o=.mod)
ZSRCS = $(ZOBJS:.o=.c)

all:
	$(MAKE) modules
	@echo "Complete making $(LIBJCUP)."

modules: $(COBJS)
#	$(INSTALL) $(ZOBJS) ../
	$(INSTALL) $(COBJS) ../
	$(INSTALL) *mod ../

clean:
	$(RM) -f *.o *.mod *.lst *.L

#$(ZOBJS): $(ZSRCS)
#	$(CC) -c $<

.SUFFIXES:
.SUFFIXES: .o .f90 .mod

.f90.o:
	$(FC) $(FFLAGS) -c $<

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)

