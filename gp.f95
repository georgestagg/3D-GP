program gp
	use params
	implicit none
	character(len=80) fname
	double precision :: ret
	CALL init_params

	!Initialise
	GRID = 0
	TIME = 0.0d0
	VOB = DBLE(VOBS)/VOBSCALE
	DT = -EYE*DTSIZE
	call calc_OBJPOT
	call approx
	write(fname, '(a,i0)') 'utils.',VOBS
	open (8, FILE = fname)
	!Iterate
	call runit(ISTEPS,0,0)
	DT = DTSIZE
	call add_noise
	call runit(NSTEPS,1,1)
	close(8)

end PROGRAM gp

subroutine runit(steps,rt,plot)
	use params
	implicit none
	integer :: steps,rt,i,plot
	if (rt == 1) then
		if (plot == 1) then
			call dump_wavefunction(0) !dump initial condition
		end if
	end if
	do i = 1, steps
		call iterate(rt)
		TIME = TIME + dble(DT)
		if (modulo(i,dumputil) == 0) then
			if (plot == 1) then
				call calc_misc
			end if
		end if
		if (modulo(i,10) == 0) then
			open (10, FILE = "STATUS")
			if (rt == 1) then
					write (unit=10,fmt="(a,i3,a,f6.2,a)")&
						"Simulating: ",VOBS, ":     ",(dble(i)/dble(steps))*100.0d0,"%"
				else
					write (unit=10,fmt="(a,i3,a,f6.2,a)")&
						"Ground State: ",VOBS, ":     " ,(dble(i)/dble(steps))*100.0d0,"%"
			end if
			close(10)
		end if
		if (modulo(i,dumpwf) == 0) then
			if (rt == 1) then
				if (plot == 1) then
					call dump_wavefunction(i)
				end if
			end if
		end if
		if(potRep .eq. 1 .and. rt == 1) then
			call calc_OBJPOT
		end if
	end do
end subroutine

subroutine dump_wavefunction (II)
	use params
	implicit none
	integer :: II,i,j,k
	character(len=80) fname
	write(fname, '(i0.4,a,i0.4)') VOBS,'.dumpwf.',II/dumpwf
	open (7, FILE = fname)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			do k = -NZ/2, NZ/2
			write (unit=7,fmt="(3f10.2,3F20.10)")&
				dble(i*DSPACE),dble(j*DSPACE),dble(k*DSPACE),dble(GRID(i,j,j)),&
				aimag(GRID(i,j,j)),DBLE(OBJPOT(i,j,j))
			end do
			write (unit=7,fmt="(a)") " "
		end do
		write (unit=7,fmt="(a)") " "
	end do
	close(7)
end subroutine
