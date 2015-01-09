
FC = gfortran

FCFLAGS = -O3

LDFLAGS = -I/usr/include/ -L/usr/lib/ -lnetcdff -lnetcdf

PROGRAMS = gp


all: $(PROGRAMS)

gp: params.o output.o bitmap.o utils.o rhs.o potential.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

%.o: %.F95
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *.in

veryclean: clean
	rm -f *~ $(PROGRAMS)
