subroutine get_vortices(vort,vortarray)
        use params    
        implicit none
        integer :: i,j,k,vort,vlength,avgd,ilineintsize
        double precision :: temp,dist,ret,psi,intervort,vortcorev,densmax,objpotcutoff
        double precision, dimension(0:NX*NY,0:2) :: vtrack,vortarray
        double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely

        vort = 0
        vortarray = 0
        vtrack = 0.0d0
        vlength = 0
        
        !vortex detecting heuristics
        if(RHSType .eq. 0) then
        	intervort = 6.0d0
        	vortcorev = 1.2d0
        	ilineintsize = NINT(intervort)
        	densmax = 0.5d0
        	objpotcutoff = 5.0d0
        end if
        if(RHSType .eq. 1) then
        	intervort = (0.3d0/DSPACE)
        	vortcorev = 1.0d0
        	ilineintsize = NINT(intervort)
        	densmax = 0.004d0
        	objpotcutoff = 50.0d0
        end if       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        call calc_phase(phase)
       	call velxy(phase,velx,vely)
       	DO i = -NX/2+10, NX/2-10
        DO j = -NY/2+10, NY/2-10
                dist = 99999999.0d0
				psi = DBLE(grid(i,j)*CONJG(grid(i,j)))
			
                if (DBLE(OBJPOT(i,j))<objpotcutoff .and. DBLE(psi)<(0.2d0*densmax)&
                		.and. SQRT((velx(i,j)**2.0d0)+(vely(i,j)**2.0d0))>vortcorev) then
                        !FOUND ONE AT i, j
                        WRITE (unit=6,fmt="(a,i0,a,i0)") "Found one at ", i,"  ",j
                        if(vlength > 0) then
                                DO k = 0, vlength-1
                                        temp = SQRT((DBLE(i-vtrack(k,0))**2.0d0) + (DBLE(j-vtrack(k,1))**2.0d0))
                                        !ITS DISTANCE FROM EXISTSING POINT k IS temp
                                        CALL LINEINTVF(velx,vely,i-ilineintsize,i+ilineintsize,j-ilineintsize,j+ilineintsize,ret)
                                        if(temp < intervort .and. sign(1.0d0,ret) .eq. sign(1.0d0,vtrack(k,2))) then
                                        	!AVERAGING vtrack(k,0),vtrack(k,1) and i,j
                                  			!WRITE (unit=6,fmt="(a, 2f10.2, a , 2i5)") "Averaging ", vtrack(k,0),vtrack(k,1), " and ", i,j
                                        	vtrack(k,0) = (vtrack(k,0)+DBLE(i))/2.0d0
                                    	    vtrack(k,1) = (vtrack(k,1)+DBLE(j))/2.0d0
                                    		!WRITE (unit=6,fmt="(a, 2f10.2)") "We get ", vtrack(k,0),vtrack(k,1)
                                    	end if
                                    	
                                    	if(temp < intervort .and. sign(1.0d0,ret) .ne. sign(1.0d0,vtrack(k,2))) then
                                        	vtrack(vlength,0) = DBLE(i)
                                        	vtrack(vlength,1) = DBLE(j)
                                        	CALL LINEINTVF(velx,vely,i-ilineintsize,i+ilineintsize,j-ilineintsize,j+ilineintsize,ret)
                             				vtrack(vlength,2) = ret
                                        	vlength = vlength + 1
                                      		!WRITE (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                               			 end if

                                        if(temp < dist) then
                                        	dist = temp
                                        end if                                        
                                end do
                                if(dist > intervort) then
                                        vtrack(vlength,0) = DBLE(i)
                                        vtrack(vlength,1) = DBLE(j)
                                        CALL LINEINTVF(velx,vely,i-ilineintsize,i+ilineintsize,j-ilineintsize,j+ilineintsize,ret)
                             			vtrack(vlength,2) = ret
                                        vlength = vlength + 1
                                      	!WRITE (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                                end if
                        else
                                vtrack(vlength,0) = DBLE(i)
                                vtrack(vlength,1) = DBLE(j)
                                CALL LINEINTVF(velx,vely,i-ilineintsize,i+ilineintsize,j-ilineintsize,j+ilineintsize,ret)
                                vtrack(vlength,2) = ret
                                vlength = vlength + 1
                                !WRITE (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                        end if
                end if
        END DO
        END DO
        vort =  vlength
		vortarray = vtrack
end subroutine


SUBROUTINE LINEINTVF(fieldx,fieldy,x,ex,y,ey,ret)
use params   
	IMPLICIT NONE
	INTEGER :: t,x,ex,y,ey
	DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: fieldx,fieldy
	DOUBLE PRECISION l1,l2,l3,l4,ret
	l1=0.0d0
	l2=0.0d0
	l3=0.0d0
	l4=0.0d0
	DO t = y,ey
		l1 = l1 + dspace*fieldy(x,t)
	END DO
	DO t = x,ex
		l2 = l2 + dspace*fieldx(t,y)
	END DO
	DO t = y,ey
		l3 = l3 + dspace*fieldy(ex,t)
	END DO
	DO t = x,ex
		l4 = l4 + dspace*fieldx(t,ey)
	END DO
	ret = l2+l3-l4-l1
END SUBROUTINE

subroutine get_delta_wrap(i,j,phase,ret)
	use params
	implicit none
	integer :: i,j
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase
	double precision :: phase1, phase2, phase3, phase4,ret
	
	call unwrap(phase(i+1,j) - phase(i,j),phase1)
	call unwrap(phase(i+1,j+1) - phase(i+1,j),phase2)
	call unwrap(phase(i,j+1) - phase(i+1,j+1),phase3)
	call unwrap(phase(i,j) - phase(i,j+1),phase4)
	ret = phase1+phase2+phase3+phase4
end subroutine

recursive subroutine unwrap(theta,ret)
	use params
	implicit none
	double precision :: theta, ret
	ret = 0.0d0
	if (abs(theta) < (PI/2.0d0)) then
		ret = theta
	else
		if (theta > (PI/2.0d0)) then
			call unwrap(theta-PI,ret)
		else
			call unwrap(theta+PI,ret)
		end if
	end if
end subroutine
