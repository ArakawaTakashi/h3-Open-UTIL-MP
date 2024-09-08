TOPDIR = ../../..
include ${TOPDIR}/Mkinclude
include files.mk

ZSRCS = $(ZOBJS:.o=.c)


all:
	echo $(ZSRCS)
	$(MAKE) modules
	@echo "Complete making $(LIBJCUP)."

modules: $(ZOBJS)
	$(INSTALL) $(ZOBJS) ../

clean:
	$(RM) -f *.o *.lst *.L

$(ZOBJS): $(ZSRCS)
#$(CC) $(CFLAGS) -c $<

.SUFFIXES:
.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c $< -I$(WAITIOINCDIR)

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)

