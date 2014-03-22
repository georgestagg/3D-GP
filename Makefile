
FC = gfortran

FCFLAGS = -O3

LDFLAGS = 

PROGRAMS = gp


all: $(PROGRAMS)

gp: params.o utils.o rhs.o potential.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

%.o: %.F95
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *.in

veryclean: clean
	rm -f *~ $(PROGRAMS)
