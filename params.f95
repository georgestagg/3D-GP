module params
  !DEFAULT VALUES---------------------------------------------------------------
  !Iterations to run
  integer :: NSTEPS=10000
  integer :: ISTEPS=2000
  integer :: VSTEPS=50

  !Resolution
  integer :: NX = 64
  integer :: NY = 64
  integer :: NZ = 64
  double precision :: DSPACE = 0.2d0, DTSIZE = 0.01d0
  
  !Dump frequency - Wavefunction - Misc Utils
  integer :: dumpwf = 100, dumputil = 100

  !GPE Type - 0 Natural Units - 1 Hamonic Oscillator Units
  integer :: RHSType = 1
  double precision :: harm_osc_C = 300.0d0
  double precision :: harm_osc_mu = 10.136d0
  complex*16 :: GAMMAC = 0.0d0
  logical :: rtNorm = .false.

  !Boundary Conditions - 0 reflective - 1 periodic
  integer :: BCX = 0
  integer :: BCY = 0
  integer :: BCZ = 0

  !Noise Amplitude - 0.001d0 works well
  integer :: noiseamp = 0.000d0

  !Flow Speed in X Dir - Start
  integer :: VOBS = 0
  double precision :: VOBSCALE = 100.0
  double precision :: DVDT = 0.0d0
  double precision :: VTVTIME = 200.0d0
 
  !Potential types - -1 none - 0 object - 3 afm-img
  logical :: enablePot = .true.
  logical :: enableTrap = .true.
  integer :: potType = -1
  integer :: trapType = 0 !traptype 0: harmonic, 1: ring
  !Enable if you need to constantly recalculate the potential
  integer :: potRep = 0

  !Object properties
  double precision :: RRX=2.0d0
  double precision :: RRY=2.0d0
  double precision :: RRZ=2.0d0
  double precision :: OBJXDASH=0.0d0
  double precision :: OBJYDASH=0.0d0
  double precision :: OBJZDASH=0.0d0
  double precision :: OBJHEIGHT=0.0d0
  !Trap
  double precision :: TXDASH=0.0d0
  double precision :: TYDASH=0.0d0
  double precision :: TZDASH=0.0d0
  double precision :: TXSCALE = 1.0d0
  double precision :: TYSCALE = 1.0d0
  double precision :: TZSCALE = 1.0d0
  double precision :: TR0 = 0.0d0

  !AFM-IMAGE
  character(2048) :: afm_filename
  integer :: afmRES = 256
  double precision :: afmXScale=1000.0d0
  double precision :: afmYscale=1000.0d0
  double precision :: afmZscale=1.0d0
  double precision :: TRUNCPARAM = 1.0d0

  !pot-image
  character(2048) :: pot_filename


  !GLOBALS----------------------------------------------------------------------
  double precision :: VOB
  double precision,parameter :: PI = 4.0d0*ATAN(1.0d0)
  complex*16 :: DT,EYE = (0.0d0,1.0d0)
  double precision :: NORM,OLDNORM = 1.0d0
  complex*16, dimension(:,:,:), ALLOCATABLE :: GRID,OBJPOT
  double precision :: TIME
  double precision :: OBJXVEL = 0.0d0,OBJYVEL = 0.0d0,OBJZVEL = 0.0d0

contains
  SUBROUTINE init_params
    IMPLICIT NONE
    afm_filename = repeat(" ", 2048) !Clear memory so entire string is blank
    pot_filename = repeat(" ", 2048) !Clear memory so entire string is blank
    include 'params.in'
    ALLOCATE(GRID(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2))
    ALLOCATE(OBJPOT(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2))
   END SUBROUTINE
end module
