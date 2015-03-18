!Runge-Kutta 4th Order!

subroutine iterate (rt)
	use params

	implicit none
	integer :: rt,BC
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: k1,k2,k3,k4

	call rhs(GRID, k1,rt)
	call rhs(GRID + 0.5d0*DT*k1,k2,rt)
	call rhs(GRID + 0.5d0*DT*k2,k3,rt)
	call rhs(GRID + DT*k3,k4,rt)
	GRID = GRID + DT*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0

	if (rt == 0 .or. rtNorm) then
		if (RHSType .eq. 0) then
			call calc_norm
			GRID = GRID/sqrt(NORM)
			GRID = GRID*sqrt(DSPACE*DSPACE*DSPACE*NX*NY*NZ)
		end if
		if (RHSType .eq. 1) then
			call calc_norm
			GRID = GRID/sqrt(NORM)
		end if
	end if
	if(rt .eq. 1 .and. abs(DVDT) > 0.0d0) then
		VOB = VOB + DVDT
	end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!

!Homogeneous Dimentionless - With Ramping!

subroutine rhs (gt, kk,rt)
	use params

	implicit none
	integer :: i,j,k,BC,rt
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: gt, kk
	kk=0
	if(RHSType .eq. 0) then
		!$OMP PARALLEL DO
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				do k = -NZ/2,NZ/2
					kk(i,j,k) = 0.5d0*(6.0d0*gt(i,j,k)- gt(BC(i,0),BC(j+1,1),BC(k,2))&
								-gt(BC(i,0),BC(j-1,1),BC(k,2))-gt(BC(i+1,0),BC(j,1),BC(k,2))&
								-gt(BC(i-1,0),BC(j,1),BC(k,2))-gt(BC(i,0),BC(j,1),BC(k-1,2))&
								-gt(BC(i,0),BC(j,1),BC(k+1,2)))/(DSPACE**2.0d0)&	!laplacian
							+gt(i,j,k)*gt(i,j,k)*CONJG(gt(i,j,k))&	!Nonlinear
					 		-gt(i,j,k)&	!Chemical Potential
							+OBJPOT(i,j,k)*gt(i,j,k)&	!Obstacle potential
							+tanh(TIME/VTVTIME)*VOB*EYE*(gt(BC(i+1,0),BC(j,1),BC(k,2))&
							-gt(BC(i-1,0),BC(j,1),BC(k,2)))/(2.0d0*DSPACE)	!Moving frame
				end do
			end do
		end do
		!$OMP END PARALLEL DO
	end if
!Harmonic Dimentionless
	if(RHSType .eq. 1) then
		!$OMP PARALLEL DO
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				do k = -NZ/2,NZ/2
				kk(i,j,k) = 0.5d0*(6.0d0*gt(i,j,k)- gt(BC(i,0),BC(j+1,1),BC(k,2))&
							-gt(BC(i,0),BC(j-1,1),BC(k,2))-gt(BC(i+1,0),BC(j,1),BC(k,2))&
							-gt(BC(i-1,0),BC(j,1),BC(k,2))-gt(BC(i,0),BC(j,1),BC(k-1,2))&
							-gt(BC(i,0),BC(j,1),BC(k+1,2)))/(DSPACE**2.0d0)&	!laplacian
						+harm_osc_C*gt(i,j,k)*gt(i,j,k)*CONJG(gt(i,j,k))&	!Nonlinear
			 			+OBJPOT(i,j,k)*gt(i,j,k)&	!potential
						- harm_osc_mu*gt(i,j,k)	!Chemical Potential
				end do
			end do
		end do
		!$OMP END PARALLEL DO
	end if
	if(DBLE(GAMMAC) > 0.0d0) then
		kk=kk/(EYE-GAMMAC)	!Damping
	else
		kk = kk/EYE !No damping
	end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function BC(s,n)
	use params
	implicit none
	integer :: s,n
	BC=s
	select case (n)
    	case (0)
    		if(s.eq.NX/2+1  .and. BCX.eq.0)BC=NX/2
			if(s.eq.NX/2+1  .and. BCX.eq.1)BC=-NX/2
    		if(s.eq.-NX/2-1 .and. BCX.eq.0)BC=-NX/2
			if(s.eq.-NX/2-1 .and. BCX.eq.1)BC=NX/2
    	case (1)
    		if(s.eq.NY/2+1  .and. BCY.eq.0)BC=NY/2
			if(s.eq.NY/2+1  .and. BCY.eq.1)BC=-NY/2
    		if(s.eq.-NY/2-1 .and. BCY.eq.0)BC=-NY/2
			if(s.eq.-NY/2-1 .and. BCY.eq.1)BC=NY/2
    	case (2)
    		if(s.eq.NZ/2+1  .and. BCZ.eq.0)BC=NZ/2
			if(s.eq.NZ/2+1  .and. BCZ.eq.1)BC=-NZ/2
    		if(s.eq.-NZ/2-1 .and. BCZ.eq.0)BC=-NZ/2
			if(s.eq.-NZ/2-1 .and. BCZ.eq.1)BC=NZ/2
	end select
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
