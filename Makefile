# Warrick's attempt at a new Makefile.
# Drawn from
# http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

#FC = g77
FC = gfortran

#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O15
#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O15 -fbounds-check
#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O15 -m32
#FFLAGS = -ffixed-line-length-none -finit-local-zero -fno-automatic -O2
#FFLAGS = -extend_source -fast # -fpe3
# FFLAGS = -extend_source -g
FFLAGS = -ffixed-line-length-none -fno-automatic -g #-fno-align-commons #-ffpe-trap=underflow #-fcheck=all -fmax-stack-var-size=999999 -fstack-arrays

#ifort options
#FFLAGS = -e -fast
#FFLAGS = -extend_source -ftz -static -C -fast

ODIR=obj
SDIR=src

_OBJ = main.o compos.o difrns.o divide.o elimn8.o equns1.o equns2.o \
funcs1.o funcs2.o nucrat.o nucrat2.o pressi.o printa.o printb.o printc.o \
remesh.o neutron.o xopac.o massloss.o diffusion.o diffusion2.o \
solver.o statef.o statel.o fdirac.o consts.o opacty.o opspln.o spline.o forceheflash.o tzo_neutron.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

VPATH = $(SDIR)

$(ODIR)/%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS)

bs: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ bs
