program gp
	use params
	
	implicit none
	character(len=80) fname
	double precision :: ret

	do LOOPNO = VOBS,VOBE,VOBST
			!Initialise
			GRID = 0
			TIME = 0.0d0
			VOB = DBLE(LOOPNO)/100.0d0
			DT = (0,-0.01d0)
			call calc_OBJPOT
			call approx
			write(fname, '(a,i0)') 'utils.',LOOPNO
			open (8, FILE = fname)
			!!!!!!!!!!!!
			call runit(2000,0,0)
			DT = (0.01d0,0)
			call add_noise
			call runit(NSTEPS,1,1)
			close(8)
	end do
end PROGRAM gp

subroutine runit(steps,rt,plot)	
	use params
	implicit none
	integer :: steps,rt,i,plot
	do i = 1, steps
		call iterate(rt)
		TIME = TIME + dble(DT)
		if(potType .eq. 2) then
			call calc_OBJPOT_osc
		end if
		if (modulo(i,10) == 0) then
			if (plot == 1) then
				call calc_misc
			end if
		end if
		if (modulo(i,100) == 0) then
			open (10, FILE = "STATUS")
			if (rt == 1) then
					write (unit=10,fmt="(a,i3,a,f6.2,a)") "Simulating: ",LOOPNO, ":     " ,(dble(i)/dble(steps))*100.0d0,"%"
				else
					write (unit=10,fmt="(a,i3,a,f6.2,a)") "Smoothing Steady Soln: ",LOOPNO, ":     " ,(dble(i)/dble(steps))*100.0d0,"%"
			end if
			close(10)
		end if
		if (modulo(i,dumpv*100) == 0) then
			if (rt == 1) then
				if (plot == 1) then	
					call dump_vortex_locations (i)
				end if	
			end if
		end if
		if (modulo(i,dumpd*100) == 0) then
			if (rt == 1) then
				if (plot == 1) then	
					call dump_density(i)
				end if	
			end if
		end if		
		if (modulo(i,dumpwf*100) == 0) then
			if (rt == 1) then
				if (plot == 1) then	
					call dump_wavefunction(i)
				end if	
			end if
		end if		
		
	end do
end subroutine



subroutine dump_density (II)
	use params
	implicit none
	integer :: II,i,j
	character(len=80) fname
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
	call calc_phase(phase)
	call velxy(phase,velx,vely)
	write(fname, '(i0,a,i0.4)') LOOPNO,'plot.',II/100
	open (7, FILE = fname)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			write (unit=7,fmt="(f10.2,f10.2,5f15.10)")&
			dble(i*DSPACE),dble(j*DSPACE),dble(GRID(i,j)*CONJG(GRID(i,j))),&
			 phase(i,j),SQRT((velx(i,j)**2.0d0)+(vely(i,j)**2.0d0)),velx(i,j),vely(i,j)
		end do
		write (unit=7,fmt="(a)") " "
	end do
	close(7)
end subroutine

subroutine dump_wavefunction (II)
	use params
	implicit none
	integer :: II,i,j
	character(len=80) fname
	write(fname, '(i0,a,i0.4)') LOOPNO,'plotwf.',II/100
	open (7, FILE = fname)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			write (unit=7,fmt="(f10.2,f10.2,5f15.10)")&
			dble(i*DSPACE),dble(j*DSPACE),dble(GRID(i,j)),aimag(GRID(i,j))
		end do
		write (unit=7,fmt="(a)") " "
	end do
	close(7)
end subroutine

subroutine dump_vortex_locations (II)
	use params
	
	implicit none
	integer :: II,k,vort
	character(len=80) fname
	double precision, dimension(0:NX*NY,0:2) :: vortarray
	vort=0
	vortarray=0.0d0
	call get_vortices(vort,vortarray)	
	write(fname, '(i0,a,i0.4)') LOOPNO,'plotvort.',II/100
	open (7, FILE = fname)
	do k=0,vort-1
		if(vortarray(k,2) > 0) then
			write (unit=7,fmt="(f10.2,f10.2,i5)") vortarray(k,0)*DSPACE,vortarray(k,1)*DSPACE,0
		else
			write (unit=7,fmt="(f10.2,f10.2,i5)"),vortarray(k,0)*DSPACE,vortarray(k,1)*DSPACE,1
		end if
	end do
	close(7)
end subroutine
