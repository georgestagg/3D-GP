!Runge-Kutta 4th Order!

subroutine iterate (rt)
	use params

	implicit none
	integer :: rt
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: k1,k2,k3,k4
	
	call rhs(GRID, k1)
	call rhs(GRID + 0.5d0*DT*k1,k2)
	call rhs(GRID + 0.5d0*DT*k2,k3)
	call rhs(GRID + DT*k3,k4)
	GRID = GRID + DT*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
	
	if (rt == 0) then	
		if (RHSType .eq. 0) then
			call calc_norm
			WRITE(6,*) NORM
			GRID = GRID/sqrt(NORM)
			GRID = GRID*sqrt(DSPACE*DSPACE*NX*NY)
		end if
		if (RHSType .eq. 1) then
			call calc_norm
			GRID = GRID/sqrt(NORM)
		end if
	end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!

!Homogeneous Dimentionless - With Ramping!

subroutine rhs (gt, kk)
	use params

	implicit none
	integer :: i,j,BC
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: gt, kk
	kk=0
	if(RHSType .eq. 0) then
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				kk(i,j) = 0.5d0*(4.0d0*gt(BC(i,0),BC(j,1))- gt(BC(i,0),BC(j+1,1))- gt(BC(i,0),BC(j-1,1))-&
						 gt(BC(i+1,0),BC(j,1))-gt(BC(i-1,0),BC(j,1)))/(DSPACE**2.0d0) + gt(i,j)*gt(i,j)*CONJG(gt(i,j))-&
						 gt(i,j) + OBJPOT(i,j)*gt(i,j)+ tanh(TIME/100.0d0)*&
						 VOB*EYE*(gt(BC(i+1,0),BC(j,1))-gt(BC(i-1,0),BC(j,1)))/(2.0d0*DSPACE)
			end do
		end do	
	end if
	
	if(RHSType .eq. 1) then
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				kk(i,j) = 0.5d0*(4.0d0*gt(BC(i,0),BC(j,1))- gt(BC(i,0),BC(j+1,1))- gt(BC(i,0),BC(j-1,1))-&
						 gt(BC(i+1,0),BC(j,1))-gt(BC(i-1,0),BC(j,1)))/(DSPACE**2.0d0) + harm_osc_C*gt(i,j)*gt(i,j)*CONJG(gt(i,j))-&
						 gt(i,j) + OBJPOT(i,j)*gt(i,j) -harm_osc_mu*gt(i,j)
			end do
		end do	
	end if	
	kk=kk/(EYE-GAMMAC)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function BC(s,n)
	use params
	implicit none
	integer :: s,n
	BC=s
	select case (n)
    	case (0)
    		if(s.eq.NX/2+1)BC=-NX/2
    		if(s.eq.-NX/2-1)BC=NX/2
    	case (1)
    		if(s.eq.NY/2+1)BC=NY/2
    		if(s.eq.-NY/2-1)BC=-NY/2
	end select
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
