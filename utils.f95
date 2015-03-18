subroutine init_RNG
	use params
	integer :: i, un, n, irand, istat
	integer, allocatable :: seed(:)
	call random_seed(size = n)
	allocate(seed(n))
	open(newunit=un, file="/dev/urandom", access="stream", &
	form="unformatted", action="read", status="old", iostat=istat)
	if (istat == 0) then
		read(un) seed
		close(un)
	else
		write (8, *) 'No /dev/urandom. Falling back to time as random seed'
		do i = 1, n
			call system_clock(seed(i))
		end do
	end if
	write (8, *) 'Initialised RNG with seed: ', seed
	call random_seed(put=seed)
	flush(8)
end subroutine

subroutine add_noise
	use params
	implicit none
	double precision :: rn1,rn2
	integer :: i,j,k,seed
	call random_number(rn1)
	call random_number(rn2)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			do k = -NZ/2, NZ/2
			GRID(i,j,k) = GRID(i,j,k) + GRID(i,j,k)*CMPLX((rn1*noiseamp)&
									-(noiseamp/2.0d0),(rn2*noiseamp)-(noiseamp/2.0d0))
			end do
		end do
	end do
end subroutine


!Approx - homogeneous!
subroutine initCond
	use params
	implicit none
	integer :: i,j,k
	if(initialCondType .eq. 0) then
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
	else if(initialCondType .eq. 1) then
		call makeRandomPhase
	end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!

subroutine makeRandomPhase
	use params
	implicit none
	double precision :: phi, rand
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: a
	integer :: i, j, k, ii2, jj2, kk2

	call calc_rpKC_rpAMP

	write (8, *) 'Complex amplitude=', rpAMP
	write (8, *) 'Cutoff wavenumber=', rpKC
	flush(8)

	do k=-NZ/2,NZ/2
		if (k <= 0) then
			kk2 = k**2
		else
			kk2 = (NZ/2-k+1)**2
		end if
		do j=-NY/2,NY/2
			if (j <= 0) then
				jj2 = j**2
			else
				jj2 = (NY/2-j+1)**2
			end if
			do i=-NX/2,NX/2
				if (i <= 0) then
					ii2 = i**2
				else
					ii2 = (NX/2-i+1)**2
				end if
				if ((ii2 + jj2 + kk2) <= rpKC) then
					call random_number(phi)
					a(i,j,k) = rpAMP*exp(2.0*PI*EYE*phi)
				else
					a(i,j,k) = 0.0d0
				end if
			end do
		end do
	end do

	call fft(a, GRID, 'backward')

end subroutine

subroutine calc_rpKC_rpAMP
	use params
	implicit none
	double precision :: ev
	ev = ((0.5*NX/PI)**2.0)*ENERV
	rpAMP = sqrt(((0.11095279993106999d0*NV**2.5d0)*(NX*NY*NZ)) / (ev**1.5d0))
	rpKC = (1.666666666666666d0*ev)/NV
end subroutine



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
	!$OMP PARALLEL DO
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
	!$OMP END PARALLEL DO
	energy = energy*DSPACE*DSPACE*DSPACE
end subroutine

subroutine insert_vortex_line(xloc,yloc,zloc,circ,rot,amp1,amp2,kk,ll,k3,dtheta)
	use params
	implicit none
	integer :: i,j,k,rot,circ
	double precision :: xloc,yloc,zloc,amp1,amp2,kk,ll,dtheta,rs,xx,yy,zz,cd,sd,k3
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: R
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: phse

	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			do k = -NZ/2, NZ/2
				if (rot .eq. 0) then
					xx = (i*DSPACE) - xloc
					yy = (j*DSPACE) - yloc
					zz = (k*DSPACE) - zloc
				else if (rot .eq. 1) then
					yy = (i*DSPACE) - xloc
					zz = (j*DSPACE) - yloc
					xx = (k*DSPACE) - zloc
				else if (rot .eq. 2) then
					zz = (i*DSPACE) - xloc
					xx = (j*DSPACE) - yloc
					yy = (k*DSPACE) - zloc
				end if
				if(abs(zz) < k3) then
					cd = amp1*cos(2.0d0*PI*zz/kk)
					sd = amp2*sin(2.0d0*PI*zz/ll)
				else
					cd = 0
					sd = 0
				end if

				rs=(xx-sd)**2.0d0 + (yy-cd)**2.0d0;
				phse(i,j,k) = exp(circ*EYE*(ATAN2(yy-cd,xx-sd) + dtheta*zz))
				R(i,j,k) = sqrt(rs*(0.3437+0.0286*rs)/(1+0.3333*rs+0.0286*rs*rs));
			end do
		end do
	end do
	GRID = GRID*R*phse;
end subroutine

subroutine insert_vortex_ring(x0,y0,z0,r0,dir,rot)
	use params
	implicit none
	integer :: i,j,k,dir,rot,n,m
	double precision :: x0,y0,z0,r0,xx,yy,zz
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: vortex_ring
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2) :: rr1,rr2,d1,d2
	double precision, dimension(:,:), ALLOCATABLE :: s

	if (rot .eq. 0) then
		ALLOCATE(s(-NX/2:NX/2,-NY/2:NY/2))
		do j=-NY/2,NY/2
		  do i=-NX/2,NX/2
		  	xx = (i*DSPACE) - x0
			yy = (j*DSPACE) - y0
			s(i,j) = sqrt(xx**2 + yy**2)
		  end do
		end do
	else if (rot .eq. 1) then
		ALLOCATE(s(-NZ/2:NZ/2,-NX/2:NX/2))
		do i=-NX/2,NX/2
		  do k=-NZ/2,NZ/2
		  	yy = (i*DSPACE) - x0
			xx = (k*DSPACE) - z0
			s(k,i) = sqrt(xx**2 + yy**2)
		  end do
		end do
	else if (rot .eq. 2) then
		ALLOCATE(s(-NY/2:NY/2,-NZ/2:NZ/2))
		do k=-NZ/2,NZ/2
		  do j=-NY/2,NY/2
			xx = (j*DSPACE) - y0
			yy = (k*DSPACE) - z0
			s(j,k) = sqrt(xx**2 + yy**2)
		  end do
		end do
	end if

    do k=-NZ/2,NZ/2
      do j=-NY/2,NY/2
        do i=-NX/2,NX/2
          if (rot .eq. 0) then
          	n=i
          	m=j
			zz = (k*DSPACE) - z0
          else if (rot .eq. 1) then
            n=k
          	m=i
			zz = (j*DSPACE) - y0
          else if (rot .eq. 2) then
            n=j
            m=k
			zz = (i*DSPACE) - x0
          end if
          d1(i,j,k) = sqrt( (zz)**2 + (s(n,m)+r0)**2 )
          d2(i,j,k) = sqrt( (zz)**2 + (s(n,m)-r0)**2 )
        end do
      end do
    end do

    rr1 = sqrt( ((0.3437d0+0.0286d0*d1**2)) / &
      (1.0d0+(0.3333d0*d1**2)+(0.0286d0*d1**4)) )
    rr2 = sqrt( ((0.3437d0+0.0286d0*d2**2)) / &
      (1.0d0+(0.3333d0*d2**2)+(0.0286d0*d2**4)) )

    do k=-NZ/2,NZ/2
      do j=-NY/2,NY/2
        do i=-NX/2,NX/2
          if (rot .eq. 0) then
          	n=i
          	m=j
			zz = (k*DSPACE) - z0
          else if (rot .eq. 1) then
            n=k
          	m=i
			zz = (j*DSPACE) - y0
          else if (rot .eq. 2) then
            n=j
            m=k
			zz = (i*DSPACE) - x0
          end if
          vortex_ring(i,j,k) =  rr1(i,j,k)*(zz+dir*EYE*(s(n,m)+r0)) * &
                                rr2(i,j,k)*(zz-dir*EYE*(s(n,m)-r0))
        end do
      end do
    end do
	GRID = GRID*vortex_ring
end subroutine

subroutine fft(in_var, out_var, dir)
	! Calculate the FFT (or inverse) of a variable
	use params
	use omp_lib
	use, intrinsic :: iso_c_binding
	implicit none
	include 'fftw3.f03'


	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2), intent(inout) :: in_var
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2,-NZ/2:NZ/2), intent(out) :: out_var
	character(*), intent(in) :: dir
	integer :: i, j, k
	integer iret
	integer*8 :: plan
	call dfftw_init_threads(iret)
	call dfftw_plan_with_nthreads(omp_get_max_threads())

	select case (dir)
	case ('forward')
		call dfftw_plan_dft_3d(plan, NX+1, NY+1, NZ+1, in_var, out_var, FFTW_FORWARD, FFTW_ESTIMATE)
	case ('backward')
		call dfftw_plan_dft_3d(plan, NX+1, NY+1, NZ+1, in_var, out_var, FFTW_BACKWARD, FFTW_ESTIMATE)
	case default
		write(6,*) 'ERROR: Not a valid direction for FFT!'
	end select

	call dfftw_execute_dft(plan, in_var, out_var)
	call dfftw_destroy_plan(plan)

	! Renormalise
	out_var = out_var / sqrt(dble((NX+1)*(NY+1)*(NZ+1)))
end subroutine fft
