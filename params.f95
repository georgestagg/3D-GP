module params
integer,parameter :: NX = 128,NY=128
double precision :: DSPACE = 0.2d0,VOB,EPS=3.0d0,R1=0.5d0
integer :: VOBS = 80, VOBE = 80, VOBST = 4,LOOPNO,OBJXDASH=0,NSTEPS=30000
double precision :: PI = 4.0d0*ATAN(1.0d0)
complex*16 :: DT,EYE = (0.0d0,1.0d0)

double precision :: NORM
complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: GRID,OBJPOT
double precision :: TIME
double precision, dimension(2) :: FVECOLD = 0.0d0
end module
