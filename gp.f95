program gp
	use params
	
	implicit none
	character(len=80) fname
	
	do LOOPNO = VOBS,VOBE,VOBST
			!Initialise
			GRID = 0
			TIME = 0.0d0
			DT = (0,-0.01d0)
			!!!!!!!!!!!!!

			write(fname, '(a,i0)') 'utils.',LOOPNO
			open (8, FILE = fname)
			call calc_OBJPOT(1)
			call approx
			call runit(2000,0,0)
			VOB=LOOPNO/100.0d0
			DT = (0.01d0,0)
			call add_noise
			call runit(NSTEPS,1,1)
	end do
end PROGRAM gp

subroutine runit(steps,rt,plot)	
	use params
	
	implicit none
	integer :: steps,rt,i,plot
	do i = 1, steps
		call iterate(rt)
		TIME = TIME + dble(DT)
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
			if (rt == 1) then
				if (plot == 1) then	
				call dump_vortex_locations (i)
				end if	
			end if
		end if
		if (modulo(i,1000) == 0) then
			call dump_density(i)
			!call dump_wavefunction(i)
		end if
	end do
end subroutine



subroutine dump_density (II)
	use params
	
	implicit none
	doUBLE PRECISION :: small
	integer :: II,i,j,vort=0
	character(len=80) fname
	complex*16 :: u,uux,uuy,uu,cuux,cuuy
	doUBLE PRECISION, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
	FORALL(i = -NX/2:NX/2, j = -NY/2:NY/2)
	phase(i,j) = ATAN(AIMAG(GRID(i,j))/dble(GRID(i,j)+TINY(0.0d0)))
	end FORALL
	call velxy(phase,velx,vely)
	write(fname, '(i0,a,i0.4)') LOOPNO,'plot.',II/100
	open (7, FILE = fname)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			write (unit=7,fmt="(f10.2,f10.2,5f15.10)")&
			dble(i),dble(j),dble(GRID(i,j)*CONJG(GRID(i,j))),&
			 phase(i,j),SQRT((velx(i,j)**2.0d0)+(vely(i,j)**2.0d0)),velx(i,j),vely(i,j)
		end do
		write (unit=7,fmt="(a)") " "
	end do
	close(7)

end subroutine

subroutine dump_vortex_locations (II)
	use params
	
	implicit none
	doUBLE PRECISION :: ret
	integer :: II,ti,tj,i,j,k,vort
	character(len=80) fname
	complex*16 :: u,uux,uuy,uu,cuux,cuuy
	doUBLE PRECISION, dimension(-NX/2:NX/2,-NY/2:NY/2) :: vorticity
	doUBLE PRECISION, dimension(0:NX*NY,0:2) :: vortarray
	vort=0
	vortarray=0.0d0
	vorticity=0.0d0
	!call GETVORTICES(vort,vortarray)	
	write(fname, '(i0,a,i0.4)') LOOPNO,'plotvort.',II/100
	open (7, FILE = fname)
	do k=0,vort-1
	if(vortarray(k,2) > 0) then
		write (unit=7,fmt="(f10.2,f10.2,i5)") vortarray(k,0),vortarray(k,1),0
	else
		write (unit=7,fmt="(f10.2,f10.2,i5)"),vortarray(k,0),vortarray(k,1),1
	end if
	end do
	close(7)

end subroutine


subroutine LINEINTVF(fieldx,fieldy,x,ex,y,ey,ret)
use params
	implicit none
	integer :: t,x,ex,y,ey
	doUBLE PRECISION, dimension(-NX/2:NX/2,-NY/2:NY/2) :: fieldx,fieldy
	doUBLE PRECISION l1,l2,l3,l4,ret
	l1=0.0d0
	l2=0.0d0
	l3=0.0d0
	l4=0.0d0
	do t = y,ey
		l1 = l1 + DSPACE*fieldy(x,t)
	end do
	do t = x,ex
		l2 = l2 + DSPACE*fieldx(t,y)
	end do
	do t = y,ey
		l3 = l3 + DSPACE*fieldy(ex,t)
	end do
	do t = x,ex
		l4 = l4 + DSPACE*fieldx(t,ey)
	end do
	ret = l2+l3-l4-l1
end subroutine

subroutine velxy(phase,velx,vely)
use params
	implicit none
	integer :: i,j
	doUBLE PRECISION, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
	doUBLE PRECISION :: temp1,temp2
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

subroutine GETVORTICES(vort,vortarray)
        use params
        
        implicit none
        integer :: i,j,k,vort,dvort,vlength,np,nn
        doUBLE PRECISION, dimension(-NX/2:NX/2,-NY/2:NY/2) :: vela
        complex*16 :: uux,uuy,uu,eng
        doUBLE PRECISION :: dist,temp,psi,ret
        doUBLE PRECISION, dimension(0:NX*NY,0:2) :: vtrack,vortarray
        doUBLE PRECISION, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
        vtrack = 0.0d0
        vlength = 0
        np=0
        nn=0
		FORALL(i = -NX/2:NX/2, j = -NY/2:NY/2)
			phase(i,j) = ATAN(AIMAG(GRID(i,j))/dble(GRID(i,j)+TINY(0.0d0)))
		end FORALL
	
		call velxy(phase,velx,vely)
        do i = -NX/2+10, NX/2-10
        do j = -NY/2+10, NY/2-10
                dist = 99999999.0d0
				uu=GRID(i,j)
				uux=(GRID(i+1,j)-GRID(i-1,j))/(2.0d0*DSPACE)
				uuy=(GRID(i,j+1)-GRID(i,j-1))/(2.0d0*DSPACE)
                eng = 0.5d0*(uux+uuy)*conjg(uux+uuy) &
				+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
				+ OBJPOT(i,j)*uu*conjg(uu)
				psi = dble(GRID(i,j)*CONJG(GRID(i,j)))
			
                if (dble(OBJPOT(i,j))<1.0d0 .and. dble(psi)<0.2d0 .and. SQRT((velx(i,j)**2.0d0)+(vely(i,j)**2.0d0))>1.2d0) then
                        !FOUND ONE AT i, j
                        if(vlength > 0) then
                                do k = 0, vlength-1
                                        temp = SQRT((dble(i-vtrack(k,0))**2.0d0) + (dble(j-vtrack(k,1))**2.0d0))
                                        !ITS DISTANCE FROM EXISTSING POINT k IS temp
                                        if(temp < 6.0d0) then
                                        !AVERAGING vtrack(k,0),vtrack(k,1) and i,j
                                  		!write (unit=6,fmt="(a, 2f10.2, a , 2i5)") "Averaging ", vtrack(k,0),vtrack(k,1), " and ", i,j
                                        	vtrack(k,0) = (vtrack(k,0)+dble(i))/2.0d0
                                    	    vtrack(k,1) = (vtrack(k,1)+dble(j))/2.0d0
                                    	!write (unit=6,fmt="(a, 2f10.2)") "We get ", vtrack(k,0),vtrack(k,1)
                                    	end if
                                        if(temp < dist) then
                                        dist = temp
                                        end if
                                end do
                                if(dist > 6.0d0 ) then
                                        vtrack(vlength,0) = dble(i)
                                        vtrack(vlength,1) = dble(j)
                                        call LINEINTVF(velx,vely,i-5,i+5,j-5,j+5,ret)
                             			vtrack(vlength,2) = ret
                                        vlength = vlength + 1
                                      	!write (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                                end if
                        else
                                vtrack(vlength,0) = dble(i)
                                vtrack(vlength,1) = dble(j)
                                call LINEINTVF(velx,vely,i-5,i+5,j-5,j+5,ret)
                                vtrack(vlength,2) = ret
                                vlength = vlength + 1
                                !write (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                        end if
                end if
        end do
        end do
        vort =  vlength
		vortarray = vtrack
        !write (unit=6,fmt="(i5)") vlength
end subroutine
