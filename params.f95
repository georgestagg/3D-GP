module params
integer,parameter :: NX = 1024,NY=512
double precision :: DSPACE = 0.4d0,VOB,RRX=6.0d0,RRY=6.0d0
integer :: VOBS = 35, VOBE = 35, VOBST = 1,LOOPNO=5,OBJXDASH=0,NSTEPS=300000
double precision :: PI = 4.0d0*ATAN(1.0d0)
complex*16 :: DT,EYE = (0.0d0,1.0d0)

double precision :: NORM
complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: GRID,OBJPOT
double precision :: TIME
double precision, dimension(2) :: FVECOLD = 0.0d0
end module
