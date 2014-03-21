module params
double precision :: VOB
integer :: LOOPNO=5
double precision,parameter :: PI = 4.0d0*ATAN(1.0d0)
complex*16 :: DT,EYE = (0.0d0,1.0d0)
include 'params.in'
double precision :: NORM
complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: GRID,OBJPOT
double precision :: TIME
double precision, dimension(2) :: FVECOLD = 0.0d0
double precision :: OBJYVEL = 0.0d0,OBJXVEL = 0.0d0
end module
