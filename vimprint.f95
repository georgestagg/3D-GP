subroutine vortex_line(circ,vlen,thetax,thetay,thetaz,xloc,yloc,zloc,wff)
	use params
	implicit none
	integer :: i,j,k,circ,BN
	double precision :: vlen,a,b,c,thetax,thetay,thetaz,phi,xloc,yloc,zloc
	complex*16, dimension(:,:,:), allocatable :: wfv
	
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: wff
	double precision, dimension(1:3,1:3) :: rx,ry,rz
	double precision, dimension(1:3) :: xyz,abc
	
	
	BN = ceiling(sqrt(DBLE(NX*NX + NY*NY + NZ*NZ)))
	ALLOCATE(wfv(-BN:BN,-BN:BN,-BN:BN))
	
	do i = -BN, BN
		do j = -BN, BN
			do k = -BN, BN
				xyz(1) = (i*DSPACE)
				xyz(2) = (j*DSPACE)
				xyz(3) = (k*DSPACE)
				if(abs(xyz(3)) < vlen/2.0d0) then
					wfv(i,j,k) = exp(circ*EYE*(ATAN2(xyz(2),xyz(1))))
				else
					wfv(i,j,k) = 1.0d0				
				end if
			end do
		end do
	end do
	
	
	rx(1,1) = 1.0d0
	rx(2,1) = 0.0d0
	rx(3,1) = 0.0d0
	rx(1,2) = 0.0d0
	rx(2,2) = cos(thetax)
	rx(3,2) = sin(thetax)
	rx(1,3) = 0.0d0
	rx(2,3) = -sin(thetax)
	rx(3,3) = cos(thetax)
	
	ry(1,1) = cos(thetay)
	ry(2,1) = 0.0d0
	ry(3,1) = -sin(thetay)
	ry(1,2) = 0.0d0
	ry(2,2) = 1.0d0
	ry(3,2) = 0.0d0
	ry(1,3) = sin(thetay)
	ry(2,3) = 0.0d0
	ry(3,3) = cos(thetay)
	
	rz(1,1) = cos(thetaz)
	rz(2,1) = sin(thetaz)
	rz(3,1) = 0.0d0
	rz(1,2) = -sin(thetaz)
	rz(2,2) = cos(thetaz)
	rz(3,2) = 0.0d0
	rz(1,3) = 0.0d0
	rz(2,3) = 0.0d0
	rz(3,3) = 1.0d0
	!rotate and translate to make final wavefunction
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			do k = -NZ/2, NZ/2
				xyz(1) = (i*DSPACE)-xloc
				xyz(2) = (j*DSPACE)-yloc
				xyz(3) = (k*DSPACE)-zloc
				!rotate
				abc = MATMUL(rz,MATMUL(ry,MATMUL(rx,xyz)))
				!translate
				wff(i,j,k) = wff(i,j,k)*wfv(nint(abc(1)/DSPACE),nint(abc(2)/DSPACE),nint(abc(3)/DSPACE))
								
			end do
		end do
	end do
	
end subroutine
