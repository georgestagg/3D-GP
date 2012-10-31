
FC = gfortran

FCFLAGS = -O3

LDFLAGS = 

PROGRAMS = gp


all: $(PROGRAMS)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

%.o: %.F95
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
