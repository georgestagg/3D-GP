subroutine add_noise
	use params
	implicit none
	integer :: i,j,k,seed
	call SYSTEM_CLOCK(COUNT=seed)
	call srand(seed)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			do k = -NZ/2, NZ/2
			GRID(i,j,k) = GRID(i,j,k) + GRID(i,j,k)*CMPLX((RAND()*noiseamp)&
									-(noiseamp/2.0d0),(RAND()*noiseamp)-(noiseamp/2.0d0))
			end do
		end do
	end do
end subroutine


!Approx - homogeneous!
subroutine approx
	use params
	implicit none
	integer :: i,j,k
	if (RHSType .eq. 0) then
		GRID =  1.00d0
	end if
	if (RHSType .eq. 1) then
		do i = -NX/2, NX/2
			do j = -NY/2, NY/2
				do k = -NZ/2, NZ/2
					GRID(i,j,k) = (harm_osc_mu - OBJPOT(i,j,k))/harm_osc_C
					if(DBLE(GRID(i,j,k)) < 0) GRID(i,j,k) = 0.0d0
				end do
			end do
		end do
		GRID = SQRT(GRID)
	end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_norm
	use params
	implicit none
	integer :: i,j
	double precision :: int_grid_3D
	NORM=int_grid_3D(DBLE(GRID*CONJG(GRID)))
	!write(6,*)NORM
end subroutine

function int_grid_3D(intThing)
	use params
	implicit none
	integer :: k,l
	double precision ::  simpsons_int_grid,int_grid_3D,a(3),b
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: intThing

	a(1) = simpsons_int_grid(intThing(:,:,-NZ/2))
	a(2) = 0.0d0
	do k = -NZ/2+1, NZ/2-1
		l = k + NZ/2
		if(modulo(l,2) == 0) then
			b = 2
		else 
			b = 4
		end if
		a(2) = a(2) + b*simpsons_int_grid(intThing(:,:,k))
	end do
	a(3) = simpsons_int_grid(intThing(:,:,-NZ/2))
	int_grid_3D = (DSPACE/3.0d0)*(a(1)+a(2)+a(3))
end function

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
	call calc_energy(energy)
    write (unit=8,fmt="(f7.2,f15.8)") time,energy
    flush(8)
end subroutine

subroutine calc_energy(energy)
	use params
	implicit none
	integer :: i,j,k
	COMPLEX*16 :: uux,uuy,uuz,uu
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
			do k = -NZ/2+1, NZ/2-1
				uu=GRID(i,j,k)
				uux=(GRID(i+1,j,k)-GRID(i-1,j,k))/(2.0d0*DSPACE)
				uuy=(GRID(i,j+1,k)-GRID(i,j-1,k))/(2.0d0*DSPACE)
				uuz=(GRID(i,j,k+1)-GRID(i,j,k-1))/(2.0d0*DSPACE)

				energy = energy + 0.5d0*((uux*CONJG(uux))+(uuy*CONJG(uuy))+(uuz*CONJG(uuz))) &
					+ 0.5d0*gg*uu*conjg(uu)*uu*conjg(uu)&
					+ OBJPOT(i,j,k)*uu*conjg(uu)
			end do
		end do
	end do
	energy = energy*DSPACE*DSPACE*DSPACE
end subroutine

subroutine insert_vortex_line(xloc,yloc,zloc,rot,circ,dtheta)
	use params
	implicit none
	integer :: i,j,k,rot,circ
	double precision :: xloc,yloc,zloc,dtheta,rs
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: R
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: phse

	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			do k = -NZ/2, NZ/2
				if (rot .eq. 0) then
					rs=(i*DSPACE)*(i*DSPACE) + (j*DSPACE)*(j*DSPACE);
					phse(i,j,k) = exp(EYE*(ATAN2((j*DSPACE),(i*DSPACE)) + dtheta*(k*DSPACE)))
				end if
				if (rot .eq. 1) then
					rs=(i*DSPACE)*(i*DSPACE) + (k*DSPACE)*(k*DSPACE)
					phse(i,j,k) = exp(EYE*(ATAN2((k*DSPACE),(i*DSPACE)) + dtheta*(j*DSPACE)))
				end if
				if (rot .eq. 2) then
					rs=(j*DSPACE)*(j*DSPACE) + (k*DSPACE)*(k*DSPACE)
					phse(i,j,k) = exp(EYE*(ATAN2((k*DSPACE),(j*DSPACE)) + dtheta*(i*DSPACE)))
				end if
				R(i,j,k) = sqrt(rs*(0.3437+0.0286*rs)/(1+0.3333*rs+0.0286*rs*rs)); 
			end do
		end do
	end do
	GRID = GRID*R*phse;
end subroutine

subroutine insert_vortex_ring(xloc,yloc,zloc,rad,dir,amp1,amp2,kk,ll,dtheta)
	use params
	implicit none
	integer :: i,j,k
	double precision :: xloc,yloc,zloc,rad,dir,amp1,amp2,kk,ll,dtheta
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: vortex_ring
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: rr1,rr2,d1,d2,pp
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: s,theta,dist
	
	do j=-NY/2, NY/2
		do i=-NX/2, NX/2
			theta(i,j) = atan2((j*DSPACE)-yloc,(i*DSPACE)-xloc)
			s(i,j) = sqrt(((i*DSPACE)-xloc)**2 + ((j*DSPACE)-yloc)**2) - amp1*sin(kk*theta(i,j))
		end do
	end do

	dist = amp2*cos(ll*theta)

	do k=-NZ/2, NZ/2
		do j=-NY/2, NY/2
			do i=-NX/2, NX/2
			pp(i,j,k) = (k*DSPACE) - amp1*cos(kk*theta(i,j))
			d1(i,j,k) = sqrt( ((pp(i,j,k)-zloc))**2 + ((s(i,j)+rad-dist(i,j)))**2 )
			d2(i,j,k) = sqrt( ((pp(i,j,k)-zloc))**2 + ((s(i,j)-rad-dist(i,j)))**2 )
			end do
		end do
	end do
	rr1 = sqrt( ((0.3437d0+0.0286d0*d1**2)) / (1.0+(0.3333d0*d1**2)+(0.0286d0*d1**4)) )
 	rr2 = sqrt( ((0.3437d0+0.0286d0*d2**2)) / (1.0+(0.3333d0*d2**2)+(0.0286d0*d2**4)) )

	do k=-NZ/2, NZ/2
		do j=-NY/2, NY/2
			do i=-NX/2, NX/2
			vortex_ring(i,j,k) = &
			rr1(i,j,k)*(((pp(i,j,k)-zloc)) + dir*EYE*((s(i,j)+rad-dist(i,j)))) * &
			rr2(i,j,k)*(((pp(i,j,k)-zloc)) - dir*EYE*((s(i,j)-rad-dist(i,j))))
			end do
		end do
	end do

	GRID = GRID*vortex_ring
end subroutine
