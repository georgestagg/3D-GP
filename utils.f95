subroutine add_noise
	use params

	implicit none
	integer :: i,j,seed
	
	call SYSTEM_CLOCK(COUNT=seed)
	call srand(seed)
	do i = -NX/2, NX/2
	do j = -NY/2, NY/2
		GRID(i,j) = GRID(i,j) + CMPLX((RAND()/100.0d0)-(1.0d0/200.0d0),(RAND()/100.0d0)-(1.0d0/200.0d0))	
	end do
	end do
end subroutine


!Approx - homogeneous!

subroutine approx		
	use params
	implicit none
	GRID =  1.00d0
end subroutine
!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_norm
	use params

	implicit none
	integer :: i,j
	NORM = 0.0d0
	do i = -NX/2,NX/2
		do j = -NY/2,NY/2
			NORM = NORM + GRID(i,j)*CONJG(GRID(i,j))		
		end do
	end do
	NORM=NORM*DSPACE*DSPACE
end subroutine

subroutine calc_misc
	use params

	implicit none
	double precision :: energy
	double precision, dimension(2) :: force

	call calc_force(force)
	call calc_energy(energy)
    write (unit=8,fmt="(f7.2,f15.8,5f15.10)") time,energy,force(1),force(2)
    flush(8)
end subroutine

subroutine calc_force(force)
	use params

	implicit none
	integer :: i,j
	COMPLEX*16 :: uux,uuy,uu
	double precision, dimension(2) :: force
	double precision, dimension(2) :: fvec
	fvec = 0.0d0

	do i = -NX/2+1, NX/2-1
		do j =  -NY/2+1, NY/2-1
			uu=GRID(i,j)
			uux=(GRID(i+1,j)-GRID(i-1,j))/(2.0d0*DSPACE)
			uuy=(GRID(i,j+1)-GRID(i,j-1))/(2.0d0*DSPACE)

			fvec(1) = fvec(1) + CONJG(uu)*EYE*uux
			fvec(2) = fvec(2) + CONJG(uu)*EYE*uuy
		end do
	end do
	fvec(1) = fvec(1)*DSPACE*DSPACE
	fvec(2) = fvec(2)*DSPACE*DSPACE

	force(1) = (fvec(1) - FVECOLD(1))/(2.0d0*DT*10.0d0)
	force(2) = (fvec(2) - FVECOLD(2))/(2.0d0*DT*10.0d0)

	FVECOLD(1) = fvec(1)
	FVECOLD(2) = fvec(2)

end subroutine

subroutine calc_energy(energy)
	use params

	implicit none
	integer :: i,j
	COMPLEX*16 :: uux,uuy,uu
	double precision :: energy
	energy = 0.0d0
	do i = -NX/2+1, NX/2-1
		do j = -NY/2+1, NY/2-1
			uu=GRID(i,j)
			uux=(GRID(i+1,j)-GRID(i-1,j))/(2.0d0*DSPACE)
			uuy=(GRID(i,j+1)-GRID(i,j-1))/(2.0d0*DSPACE)

			energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
				+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
				+ OBJPOT(i,j)*uu*conjg(uu)
		end do
	end do
	energy = energy*DSPACE*DSPACE
end subroutine

subroutine calc_phase(phase)
	use params
	implicit none
	integer :: i,j
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase
	forall(i = -NX/2:NX/2, j = -NY/2:NY/2)
		phase(i,j) = atan(aimag(GRID(i,j))/dble(GRID(i,j)+tiny(0.0d0)))
	end forall
end subroutine

subroutine velxy(phase,velx,vely)
	use params
	implicit none
	integer :: i,j
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
	double precision :: temp1,temp2
	velx = 0.0d0
	vely = 0.0d0
	do i = -NX/2+1,NX/2-1
	do j = -NY/2+1,NY/2-1
		if (phase(i+1,j)-phase(i-1,j)<-(pi/2.0d0)) then
			temp1 = phase(i+1,j)-(phase(i-1,j) - PI)
		else if (phase(i+1,j)-phase(i-1,j)>(pi/2.0d0)) then
			temp1 = phase(i+1,j)-(phase(i-1,j) + PI)
		else
			temp1 = phase(i+1,j)-phase(i-1,j)
		end if
		velx(i,j) = temp1
	end do
	end do
	
	do i = -NX/2+1,NX/2-1
	do j = -NY/2+1,NY/2-1
		if (phase(i,j+1)-phase(i,j-1)<-(pi/2.0d0)) then
			temp1 = phase(i,j+1)-(phase(i,j-1) - PI)
		else if (phase(i,j+1)-phase(i,j-1)>(pi/2.0d0)) then
			temp1 = phase(i,j+1)-(phase(i,j-1) + PI)
		else
			temp1 = phase(i,j+1)-phase(i,j-1)
		end if
		vely(i,j) = temp1
	end do
	end do
end subroutine