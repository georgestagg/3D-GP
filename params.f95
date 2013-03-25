module params
integer,parameter :: NX = 2048,NY=512
double precision :: DT, DSPACE = 0.4d0,VOB,EPS=3.0d0,R1=3.0d0
integer :: VOBS = 44, VOBE = 44, VOBST = 4,LOOPNO,OBJXDASH=980,NSTEPS=3000
double precision :: PI = 4.0d0*ATAN(1.0d0)
complex*16 :: EYE = (0.0d0,1.0d0)

double precision :: NORM
complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: GRID,OBJPOT
double precision :: TIME
double precision, dimension(2) :: FVECOLD = 0.0d0
end module
