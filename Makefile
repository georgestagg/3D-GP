
FC = gfortran

FCFLAGS = -O3 -fcheck=all -fopenmp

LDFLAGS = -L/usr/lib/ -L/usr/local/lib/ -lfftw3 -lm -lnetcdff -lnetcdf 

LD_LIBRARY_PATH = -I/usr/local/include/ -I/usr/include/

PROGRAMS = gp


all: $(PROGRAMS)

gp: params.o rhs.o output.o bitmap.o utils.o potential.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $< $(LD_LIBRARY_PATH)

%.o: %.F95
	$(FC) $(FCFLAGS) -c $< $(LD_LIBRARY_PATH)

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *.in

veryclean: clean
	rm -f *~ $(PROGRAMS)
