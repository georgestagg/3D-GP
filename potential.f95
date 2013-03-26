subroutine calc_OBJPOT (obj)				
	use params
	implicit none
	integer :: i,j,obj
	double precision :: rx,ry,r2
	r2 = R1*EPS
	do i = -NX/2,NX/2
		do j = -NY/2,NY/2
			rx = dble(i-OBJXDASH)
			ry = dble(j)
			OBJPOT(i,j) = 100.0d0*EXP(-(1.0d0/R1**2.0d0)*(rx*DSPACE)**2.0d0 - (1.0d0/r2**2.0d0)*(ry*DSPACE)**2.0d0)
		end do
	end do
end subroutine