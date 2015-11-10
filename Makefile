
FC = gfortran

FCFLAGS = -O3 -fopenmp

LDFLAGS = -L/usr/lib/ -L/usr/local/lib/ -lfftw3 -lm -lnetcdff -lnetcdf 

LDLIB = -I/state/partition1/apps/fftw/fftw-3.3.4/include/ -I/home/a8034837/libs/netcdf-fortran-4.4.2/include/ -I/usr/local/include/ -I/usr/include/

PROGRAMS = gp


all: $(PROGRAMS)

gp: params.o rhs.o output.o vimprint.o bitmap.o utils.o potential.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $< $(LDLIB)

%.o: %.F95
	$(FC) $(FCFLAGS) -c $< $(LDLIB)

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *.in

veryclean: clean
	rm -f *~ $(PROGRAMS)
