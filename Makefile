
FC = gfortran

FCFLAGS = -O3 -fcheck=all

LDFLAGS = -L/usr/local/lib/ -lfftw3_omp -lfftw3 -lm -lnetcdff -lnetcdf -lcurl -lhdf5 -lhdf5_hl -fopenmp 

LD_LIBRARY_PATH = -I/usr/local/include/

PROGRAMS = gp


all: $(PROGRAMS)

gp: params.o output.o bitmap.o utils.o rhs.o potential.o

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
