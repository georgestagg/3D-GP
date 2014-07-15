subroutine add_noise
	use params
	implicit none
	integer :: i,j,seed
	call SYSTEM_CLOCK(COUNT=seed)
	call srand(seed)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			GRID(i,j) = GRID(i,j) + GRID(i,j)*CMPLX((RAND()*noiseamp)&
									-(noiseamp/2.0d0),(RAND()*noiseamp)-(noiseamp/2.0d0))
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
	double precision :: simpsons_int_grid
	NORM=simpsons_int_grid(DBLE(GRID*CONJG(GRID)))
	write(6,*) NORM
end subroutine

function simpsons_int_grid(intThing)
	use params
	implicit none
	double precision :: simpsons_int_grid,a(16)
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: intThing

	a(1) = intThing(-NX/2,-NY/2)
	a(2) = intThing(-NX/2,NY/2)
	a(3) = intThing(NX/2,-NY/2)
	a(4) = intThing(NX/2,NY/2)

	a(5) = 4.0d0*sum(intThing(-NX/2,-NY/2+1:NY/2-1:2))
	a(6) = 2.0d0*sum(intThing(-NX/2,-NY/2+2:NY/2-1:2))
	a(7) = 4.0d0*sum(intThing(NX/2,-NY/2+1:NY/2-1:2))
	a(8) = 2.0d0*sum(intThing(NX/2,-NY/2+2:NY/2-1:2))

	a(9) = 4.0d0*sum(intThing(-NX/2+1:NX/2-1:2,-NY/2))
   a(10) = 2.0d0*sum(intThing(-NX/2+2:NX/2-1:2,-NY/2))
   a(11) = 4.0d0*sum(intThing(-NX/2+1:NX/2-1:2,NY/2))
   a(12) = 2.0d0*sum(intThing(-NX/2+2:NX/2-1:2,NY/2))

   a(13) = 16.0d0*sum(intThing(-NX/2+1:NX/2-1:2,-NY/2+1:NY/2-1:2))
   a(14) = 8.0d0* sum(intThing(-NX/2+1:NX/2-1:2,-NY/2+2:NY/2-1:2))
   a(15) = 8.0d0* sum(intThing(-NX/2+2:NX/2-1:2,-NY/2+1:NY/2-1:2))
   a(16) = 4.0d0*sum(intThing(-NX/2+2:NX/2-1:2,-NY/2+2:NY/2-1:2))

   simpsons_int_grid = (1.0d0/9.0d0)*(DSPACE*DSPACE)*sum(a)
end function

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

subroutine calc_force_2D(force2d)
  use params

  implicit none
  integer :: i,j
  COMPLEX*16 :: uu,uux,uuy
  double precision, dimension(2,-NX/2:NX/2,-NY/2:NY/2) :: force2d
  do i = -NX/2+1, NX/2-1
      do j = -NY/2+1, NY/2-1
        uu=GRID(i,j)*CONJG(GRID(i,j))
        uux=(OBJPOT(i+1,j)-OBJPOT(i-1,j))/(2.0d0*DSPACE)
        uuy=(OBJPOT(i,j+1)-OBJPOT(i,j-1))/(2.0d0*DSPACE)

  		  force2d(1,i,j) = uu*uux
		  	force2d(2,i,j) = uu*uuy
      end do
  end do
end subroutine

subroutine cross_2d(a, b,ret)
	implicit none
	double precision :: ret
    double precision, dimension(2) :: a,b
    ret = a(1)*b(2) - a(2)*b(1)
END subroutine

subroutine calc_net_torque(torque)
	use params
	implicit none
	integer :: i,j
	double precision, dimension(2,-NX/2:NX/2,-NY/2:NY/2) :: force2d
	double precision, dimension(2) :: r_rel
	double precision :: torque,crossp

	torque = 0.0d0
	call calc_force_2D(force2d)
	do i = -NX/2+1, NX/2-1
		do j = -NY/2+1, NY/2-1
			r_rel(1) = (dble(i)*DSPACE)-OBJXDASH
			r_rel(2) = (dble(j)*DSPACE)-OBJYDASH
			call cross_2d(r_rel,force2d(:,i,j),crossp)
		   torque = torque + crossp
		end do
	end do
end subroutine

subroutine calc_new_obj_angle
  use params
	implicit none
	double precision :: nettorque
	call calc_net_torque(nettorque)
	OBJANGLEV = OBJANGLEV + MOMINERTIA*(nettorque*DT)
	OBJANGLE  = OBJANGLE  + (OBJANGLEV*DT)
end subroutine

subroutine calc_OBJNEWTON
	use params
	implicit none
	double precision :: beta = 2.5d0, forced, mass
	double precision, dimension(2,-NX/2:NX/2,-NY/2:NY/2) :: force2d

	call calc_force_2D(force2d)
	mass = beta/(w0**2.0d0)

	forced = FVAL*SIN(DBLE(TIME)*WF)
	OBJXVEL = OBJXVEL + ((sum(force2d(1,:,:))*DSPACE*DSPACE&
						-beta*OBJXDASH + forced)/mass)*DBLE(DT)
	OBJXDASH = OBJXDASH + OBJXVEL*DBLE(DT)
	CALL calc_OBJPOT_osc
end subroutine


subroutine calc_energy(energy)
	use params
	implicit none
	integer :: i,j
	COMPLEX*16 :: uux,uuy,uu
	double precision :: energy,gg
	energy = 0.0d0
	if(RHSType .eq. 0) then
		gg=1.0d0
	end if
	if(RHSType .eq. 1) then
		gg=harm_osc_C
	end if
	do i = -NX/2+1, NX/2-1
		do j = -NY/2+1, NY/2-1
			uu=GRID(i,j)
			uux=(GRID(i+1,j)-GRID(i-1,j))/(2.0d0*DSPACE)
			uuy=(GRID(i,j+1)-GRID(i,j-1))/(2.0d0*DSPACE)

			energy = energy + 0.5d0*((uux*CONJG(uux))+(uuy*CONJG(uuy))) &
				+ 0.5d0*gg*uu*conjg(uu)*uu*conjg(uu)&
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
