module runparam
INTEGER :: s,ls=1,le=1
end module
module plane
INTEGER,PARAMETER :: NX = 2048,NY=512
DOUBLE PRECISION :: dspace = 0.4d0,PI = 4.0d0*ATAN(1.0d0), norm, vob
COMPLEX*16, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: grid,vtrap,boundary,koch=0.0d0
DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2,2) :: fvec
INTEGER, DIMENSION(0:NX*NY,0:1) :: ktrack
INTEGER :: klength = 0
COMPLEX*16 :: eye = (0.0d0,1.0d0)
END MODULE
module timee
DOUBLE PRECISION :: time
integer :: ticks=0,loopno,itime=0
END MODULE

PROGRAM gp
	use runparam
	use plane
	use timee
	IMPLICIT NONE
	INTEGER :: I,temp,j						
	COMPLEX*16 :: dt,done	
	CHARACTER(len=80) fname
	
	INTERFACE 
   		SUBROUTINE KOCHLOOP(s,oix,oiy,iter1,iter2,iter3,iter4, iter5,iter6)
     		INTEGER :: iter1,iter2,iter3,iter4
		integer :: oix,oiy
		double precision ::s,ns,ox,oy
   		END SUBROUTINE
	END INTERFACE
	DO loopno =44,44,4
			!Initialise
			koch = 0.0d0
			grid = 0
			ticks=0
			done = 0
			time = 0.0d0
			itime=1
			dt = (0,-0.01d0)
			!!!!!!!!!!!!!
			WRITE(fname, '(a,i0)') 'energy.',loopno
			OPEN (8, FILE = fname)
			ktrack = -8000000
			klength = 0
			boundary = 0

			CALL CALC_VTRAP(1)
			CALL APPROX
			CALL RUNIT(2000,dt,0,0)
			itime=0
			vob=loopno/100.0d0
			dt = (0.01d0,0)		
			CALL RUNIT(200,dt,1,0)
			CALL ADDNOISE
			CALL RUNIT(300000,dt,1,1)
	END DO
END PROGRAM gp



SUBROUTINE ADDNOISE
use plane
use runparam
use timee
IMPLICIT NONE
integer :: i,j
DO i = -NX/2, NX/2
DO j = -NY/2, NY/2
	grid(i,j) = grid(i,j) + CMPLX((RAND()/100.0d0)-(1.0d0/200.0d0),(RAND()/100.0d0)-(1.0d0/200.0d0))	
END DO
END DO
END SUBROUTINE

SUBROUTINE RUNIT(steps,dt,rt,plot)	
	use plane
	use runparam
	use timee
	IMPLICIT NONE
	INTEGER :: steps,rt,I,plot
	COMPLEX*16 :: dt
	COMPLEX*16, DIMENSION(4) :: statuss
	DO I = 1, steps
		CALL ITERATE(dt,rt)
		time = time + DBLE(dt)
		IF (MODULO(I,10) == 0) THEN
			IF (plot == 1) THEN
				CALL MAKEENERGY
			END IF
		END IF
		IF (MODULO(I,100) == 0) THEN
			OPEN (10, FILE = "STATUS")
			IF (rt == 1) THEN
					WRITE (unit=10,fmt="(a,i3,a,f6.2,a)") "Simulating: ",loopno, ":     " ,(DBLE(I)/DBLE(steps))*100.0d0,"%"
				ELSE
					WRITE (unit=10,fmt="(a,i3,a,f6.2,a)") "Smoothing Steady Soln: ",loopno, ":     " ,(DBLE(I)/DBLE(steps))*100.0d0,"%"
			END IF
			CLOSE(10)
			IF (rt == 1) THEN
				IF (plot == 1) THEN
				!CALL PLOT3D(I)	
				CALL PLOTVORTICITY (I)
				END IF	
			END IF
		END IF
		IF (MODULO(I,1000) == 0) THEN
			CALL PLOT3D(I)
		END IF
	END DO
END SUBROUTINE

SUBROUTINE ITERATE (dt,rt)
	use plane
	IMPLICIT NONE
	integer i,j,rt
	COMPLEX*16 :: dt
	COMPLEX*16, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: k1,k2,k3,k4
	
	CALL RHS(grid, k1)
	CALL RHS(grid + 0.5d0*dt*k1,k2)
	CALL RHS(grid + 0.5d0*dt*k2,k3)
	CALL RHS(grid + dt*k3,k4)
	grid = grid + dt*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
	
	IF (rt == 0) THEN	
		CALL CALC_NORM
		grid = sqrt(dspace*dspace*NX*NY)*grid/sqrt(norm)
	end if

END SUBROUTINE

SUBROUTINE CALC_VTRAP (obj)				
	use runparam
	use plane
	use timee
	IMPLICIT NONE
	INTEGER :: i,j,obj
	DOUBLE PRECISION :: T,v,xdash,rx,ry,d1,d2,eps
	COMPLEX*16, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: vobj
	!vtrap = koch
	eps = 3.0d0
	d1 = 3.0d0
	d2 = d1*eps
	DO i = -NX/2,NX/2
		DO j = -NY/2,NY/2
			rx = DBLE(i-980)
			ry = DBLE(j)
			vtrap(i,j) = 100.0d0*EXP(-(1.0d0/d1**2.0d0)*((rx)*dspace)**2.0d0 - (1.0d0/d2**2.0d0)*(ry*dspace)**2.0d0)
		END DO
	END DO
	ticks=ticks+1	
END SUBROUTINE

SUBROUTINE APPROX			
	use plane
	USE runparam
	IMPLICIT NONE
	grid =  1.00d0
END SUBROUTINE

SUBROUTINE CALC_NORM
	use plane
	IMPLICIT NONE
	INTEGER :: i,j
	norm = 0.0d0
	DO i = -NX/2,NX/2
		DO j = -NY/2,NY/2
			norm = norm + grid(i,j)*CONJG(grid(i,j))		
		END DO
	END DO
	norm=norm*dspace*dspace
END SUBROUTINE

SUBROUTINE RHS (gt, kk)
	use plane
	use timee
	IMPLICIT NONE
	INTEGER :: i,j,k,BC
	COMPLEX*16, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: gt, kk
	
	kk=0

	DO i = -NX/2,NX/2
		DO j = -NY/2,NY/2
			
			kk(i,j) = 0.5d0*(4.0d0*gt(BC(i,0),BC(j,1))- gt(BC(i,0),BC(j+1,1))- gt(BC(i,0),BC(j-1,1))-&
					 gt(BC(i+1,0),BC(j,1))-gt(BC(i-1,0),BC(j,1)))/(dspace**2.0d0) + gt(i,j)*gt(i,j)*CONJG(gt(i,j))-&
					 gt(i,j) + vtrap(i,j)*gt(i,j)+ tanh(time/100.0d0)*&
					 vob*eye*(gt(BC(i+1,0),BC(j,1))-gt(BC(i-1,0),BC(j,1)))/(2.0d0*dspace)
		END DO
	END DO
					
	kk=kk/(eye)

END SUBROUTINE

INTEGER FUNCTION BC(s,n)
	use plane
	IMPLICIT NONE
	INTEGER :: s,n
	BC=s
	SELECT CASE (n)
    case (0)
     if(s.eq.NX/2+1)BC=-NX/2
     if(s.eq.-NX/2-1)BC=NX/2
    case (1)
     if(s.eq.NY/2+1)BC=-NY/2
     if(s.eq.-NY/2-1)BC=NY/2
	END SELECT
END FUNCTION



SUBROUTINE MAKEENERGY
	use runparam
	use plane
	use timee
	IMPLICIT NONE
	INTEGER :: i,j,vort,np,nn
	COMPLEX*16, DIMENSION(4) :: energyf
	COMPLEX*16, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: vobj,vobj1,vobj2,vobj3,vobj4
	COMPLEX*16 :: uux,uuy,uu
	DOUBLE PRECISION :: T,v,force,ydash,energy,ret
	DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
	FORALL(i = -NX/2:NX/2, j = -NY/2:NY/2)
	phase(i,j) = ATAN(AIMAG(grid(i,j))/DBLE(grid(i,j)+TINY(0.0d0)))
	END FORALL
	
	CALL velxy(phase,velx,vely)
	CALL LINEINTVF(velx,vely,-30,-10,-10,10,ret)
	force = 0
	energy = 0
	DO i = -NX/2+1, NX/2-1
		DO j = -NY/2+1, NY/2-1
			uu=grid(i,j)
			uux=(grid(i+1,j)-grid(i-1,j))/(2.0d0*dspace)
			uuy=(grid(i,j+1)-grid(i,j-1))/(2.0d0*dspace)

			energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
				+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
				+ vtrap(i,j)*uu*conjg(uu)
		END DO
	END DO
	DO i = -NX/2+1, NX/2-1
		uu=grid(i,NY/2)
		uux=(grid(i+1,NY/2)-grid(i-1,NY/2))/(2.0d0*dspace)
		uuy=(grid(i,-NY/2)-grid(i,NY/2-1))/(2.0d0*dspace)

		energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
			+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
			+ vtrap(i,NY/2)*uu*conjg(uu)
		
		uu=grid(i,-NY/2)
		uux=(grid(i+1,-NY/2)-grid(i-1,-NY/2))/(2.0d0*dspace)
		uuy=(grid(i,-NY/2+1)-grid(i,NY/2))/(2.0d0*dspace)

		energy =energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
			+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
			+ vtrap(i,-NY/2)*uu*conjg(uu)
	END DO
	DO j = -NY/2+1, NY/2-1
		uu=grid(NX/2,j)
		uux=(grid(-NX/2,j)-grid(NX/2-1,j))/(2.0d0*dspace)
		uuy=(grid(NX/2,j+1)-grid(NX/2,j-1))/(2.0d0*dspace)

		energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
			+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
			+ vtrap(NX/2,j)*uu*conjg(uu)

		uu=grid(-NX/2,j)
		uux=(grid(-NX/2+1,j)-grid(NX/2,j))/(2.0d0*dspace)
		uuy=(grid(-NX/2,j+1)-grid(-NX/2,j-1))/(2.0d0*dspace)

		energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
			+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
			+ vtrap(-NX/2,j)*uu*conjg(uu)
	END DO
		uu=grid(NX/2,NY/2)
		uux=(grid(-NX/2,NY/2)-grid(NX/2-1,NY/2))/(2.0d0*dspace)
		uuy=(grid(NX/2,-NY/2)-grid(NX/2,NY/2-1))/(2.0d0*dspace)
		energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
			+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
			+ vtrap(NX/2,NY/2)*uu*conjg(uu)		
		uu=grid(NX/2,-NY/2)
		uux=(grid(-NX/2,-NY/2)-grid(NX/2-1,-NY/2))/(2.0d0*dspace)
		uuy=(grid(NX/2,-NY/2+1)-grid(NX/2,NY/2))/(2.0d0*dspace)
		energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
			+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
			+ vtrap(NX/2,-NY/2)*uu*conjg(uu)		
		uu=grid(-NX/2,NY/2)
		uux=(grid(-NX/2+1,NY/2)-grid(NX/2,NY/2))/(2.0d0*dspace)
		uuy=(grid(-NX/2,-NY/2)-grid(-NX/2,NY/2-1))/(2.0d0*dspace)
		energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
			+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
			+ vtrap(-NX/2,NY/2)*uu*conjg(uu)		
		uu=grid(-NX/2,-NY/2)
		uux=(grid(-NX/2+1,-NY/2)-grid(NX/2,-NY/2))/(2.0d0*dspace)
		uuy=(grid(-NX/2,-NY/2+1)-grid(-NX/2,NY/2))/(2.0d0*dspace)
		energy = energy + 0.5d0*(uux+uuy)*conjg(uux+uuy) &
			+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
			+ vtrap(-NX/2,-NY/2)*uu*conjg(uu)

	energy = energy*dspace*dspace
		
	DO i = -NX/2+1, NX/2-1
		DO j =  -NY/2+1, NY/2-1
			uux=(vtrap(i+1,j)-vtrap(i-1,j))/(2.0d0*dspace)
			uuy=(vtrap(i,j+1)-vtrap(i,j-1))/(2.0d0*dspace)
			fvec(i,j,1) = grid(i,j)*CONJG(grid(i,j))*uux
			fvec(i,j,2) = grid(i,j)*CONJG(grid(i,j))*uuy
		END DO
	END DO

	DO i = -NX/2+1, NX/2-1
		DO j =  -NY/2+1, NY/2-1
		force = force + fvec(i,j,1)
		END DO
	END DO
	
	vort=0
	force = -force*dspace*dspace
    CALL COUNTVORTICES(vort,np,nn)
    WRITE (unit=8,fmt="(f7.2,f15.5,f15.5,i5,i5,i5)") time,energy,force,vort,np,nn
    FLUSH(8)
END SUBROUTINE



SUBROUTINE PLOT3D (II)
	use runparam
	use plane
	use timee
	IMPLICIT NONE
	DOUBLE PRECISION :: small
	INTEGER :: II,i,j,vort=0
	CHARACTER(len=80) fname
	COMPLEX*16 :: u,uux,uuy,uu,cuux,cuuy
	DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
	FORALL(i = -NX/2:NX/2, j = -NY/2:NY/2)
	phase(i,j) = ATAN(AIMAG(grid(i,j))/DBLE(grid(i,j)+TINY(0.0d0)))
	END FORALL
	CALL velxy(phase,velx,vely)
	WRITE(fname, '(i0,a,i0.4)') loopno,'plot.',II/100
	OPEN (7, FILE = fname)
	DO i = -NX/2, NX/2
		DO j = -NY/2, NY/2
			WRITE (unit=7,fmt="(f10.2,f10.2,5f15.10)")&
			DBLE(i),DBLE(j),DBLE(grid(i,j)*CONJG(grid(i,j))),&
			 phase(i,j),SQRT((velx(i,j)**2.0d0)+(vely(i,j)**2.0d0)),velx(i,j),vely(i,j)
		END DO
		WRITE (unit=7,fmt="(a)") " "
	END DO
	CLOSE(7)

END SUBROUTINE

SUBROUTINE PLOTVORTICITY (II)
	use runparam
	use plane
	use timee
	IMPLICIT NONE
	DOUBLE PRECISION :: ret
	INTEGER :: II,ti,tj,i,j,k,vort
	CHARACTER(len=80) fname
	COMPLEX*16 :: u,uux,uuy,uu,cuux,cuuy
	DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: vorticity
	DOUBLE PRECISION, DIMENSION(0:NX*NY,0:2) :: vortarray
	vort=0
	vortarray=0.0d0
	vorticity=0.0d0
	CALL GETVORTICES(vort,vortarray)	
	WRITE(fname, '(i0,a,i0.4)') loopno,'plotvort.',II/100
	OPEN (7, FILE = fname)
	DO k=0,vort-1
	if(vortarray(k,2) > 0) then
		WRITE (unit=7,fmt="(f10.2,f10.2,i5)") vortarray(k,0),vortarray(k,1),0
	else
		WRITE (unit=7,fmt="(f10.2,f10.2,i5)"),vortarray(k,0),vortarray(k,1),1
	end if
	END DO
	CLOSE(7)

END SUBROUTINE


SUBROUTINE LINEINTVF(fieldx,fieldy,x,ex,y,ey,ret)
use plane
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

SUBROUTINE velxy(phase,velx,vely)
use plane
	IMPLICIT NONE
	INTEGER :: i,j
	DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
	DOUBLE PRECISION :: temp1,temp2
	velx = 0.0d0
	vely = 0.0d0
	DO i = -NX/2+1,NX/2-1
	DO j = -NY/2+1,NY/2-1
		if (phase(i+1,j)-phase(i-1,j)<-(pi/2.0d0)) then
			temp1 = phase(i+1,j)-(phase(i-1,j) - PI)
		else if (phase(i+1,j)-phase(i-1,j)>(pi/2.0d0)) then
			temp1 = phase(i+1,j)-(phase(i-1,j) + PI)
		else
			temp1 = phase(i+1,j)-phase(i-1,j)
		end if
		velx(i,j) = temp1
	END DO
	END DO
	
	DO i = -NX/2+1,NX/2-1
	DO j = -NY/2+1,NY/2-1
		if (phase(i,j+1)-phase(i,j-1)<-(pi/2.0d0)) then
			temp1 = phase(i,j+1)-(phase(i,j-1) - PI)
		else if (phase(i,j+1)-phase(i,j-1)>(pi/2.0d0)) then
			temp1 = phase(i,j+1)-(phase(i,j-1) + PI)
		else
			temp1 = phase(i,j+1)-phase(i,j-1)
		end if
		vely(i,j) = temp1
	END DO
	END DO
END SUBROUTINE

SUBROUTINE COUNTVORTICES(vort,np,nn)
        use runparam
        use plane
        use timee
        IMPLICIT NONE
        INTEGER :: i,j,k,vort,dvort,vlength,np,nn
        DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: vela
        COMPLEX*16 :: uux,uuy,uu,eng
        DOUBLE PRECISION :: dist,temp,psi,ret
        INTEGER, DIMENSION(0:NX*NY,0:2) :: vtrack
        DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
        vtrack = 0.0d0
        vlength = 0
        np=0
        nn=0
		FORALL(i = -NX/2:NX/2, j = -NY/2:NY/2)
			phase(i,j) = ATAN(AIMAG(grid(i,j))/DBLE(grid(i,j)+TINY(0.0d0)))
		END FORALL
	
		CALL velxy(phase,velx,vely)
        DO i = -NX/2+10, NX/2-10
        DO j = -NY/2+10, NY/2-10
                dist = 99999999.0d0
				uu=grid(i,j)
				uux=(grid(i+1,j)-grid(i-1,j))/(2.0d0*dspace)
				uuy=(grid(i,j+1)-grid(i,j-1))/(2.0d0*dspace)
                eng = 0.5d0*(uux+uuy)*conjg(uux+uuy) &
				+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
				+ vtrap(i,j)*uu*conjg(uu)
				psi = DBLE(grid(i,j)*CONJG(grid(i,j)))
			
                if (DBLE(vtrap(i,j))<0.1d0 .and. DBLE(psi)<0.2d0 .and. SQRT((velx(i,j)**2.0d0)+(vely(i,j)**2.0d0))>1.2d0) then
                        !WRITE (unit=6,fmt="(a,i4,i4)") "FOUND ONE AT ", i, j
                        if(vlength > 0) then
                                DO k = 0, vlength-1
                                        temp = SQRT((DBLE(i-vtrack(k,0))**2.0d0) + (DBLE(j-vtrack(k,1))**2.0d0))
                                        !WRITE (unit=6,fmt="(a,i4,i4,a,f10.7)") "ITS DISTANCE FROM EXISTSING POINT ", vtrack(k,0),vtrack(k,1), " IS ", temp
                                        if(temp < dist) then
                                        dist = temp
                                        end if
                                end do
                                if(dist > 6.0d0 ) then
                                        vtrack(vlength,0) = i
                                        vtrack(vlength,1) = j
                                        CALL LINEINTVF(velx,vely,i-5,i+5,j-5,j+5,ret)
                             			vtrack(vlength,2) = INT(ret)
                                        vlength = vlength + 1
                                       !WRITE (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                                end if
                        else
                                vtrack(vlength,0) = i
                                vtrack(vlength,1) = j
                                CALL LINEINTVF(velx,vely,i-5,i+5,j-5,j+5,ret)
                                vtrack(vlength,2) = INT(ret)
                                vlength = vlength + 1
                                !WRITE (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                        end if
                end if
        END DO
        END DO
        vort =  vlength
        DO i=0,vlength-1
        	if(vtrack(i,2) > 0) then
        	np = np + 1
        	else
        	nn = nn + 1
        	end if
        END DO
        !WRITE (unit=6,fmt="(i5)") vlength
END SUBROUTINE

SUBROUTINE GETVORTICES(vort,vortarray)
        use runparam
        use plane
        use timee
        IMPLICIT NONE
        INTEGER :: i,j,k,vort,dvort,vlength,np,nn
        DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: vela
        COMPLEX*16 :: uux,uuy,uu,eng
        DOUBLE PRECISION :: dist,temp,psi,ret
        DOUBLE PRECISION, DIMENSION(0:NX*NY,0:2) :: vtrack,vortarray
        DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
        vtrack = 0.0d0
        vlength = 0
        np=0
        nn=0
		FORALL(i = -NX/2:NX/2, j = -NY/2:NY/2)
			phase(i,j) = ATAN(AIMAG(grid(i,j))/DBLE(grid(i,j)+TINY(0.0d0)))
		END FORALL
	
		CALL velxy(phase,velx,vely)
        DO i = -NX/2+10, NX/2-10
        DO j = -NY/2+10, NY/2-10
                dist = 99999999.0d0
				uu=grid(i,j)
				uux=(grid(i+1,j)-grid(i-1,j))/(2.0d0*dspace)
				uuy=(grid(i,j+1)-grid(i,j-1))/(2.0d0*dspace)
                eng = 0.5d0*(uux+uuy)*conjg(uux+uuy) &
				+ 0.5d0*((uu*conjg(uu)-1.0d0)**2.0d0)&
				+ vtrap(i,j)*uu*conjg(uu)
				psi = DBLE(grid(i,j)*CONJG(grid(i,j)))
			
                if (DBLE(vtrap(i,j))<1.0d0 .and. DBLE(psi)<0.2d0 .and. SQRT((velx(i,j)**2.0d0)+(vely(i,j)**2.0d0))>1.2d0) then
                        !FOUND ONE AT i, j
                        if(vlength > 0) then
                                DO k = 0, vlength-1
                                        temp = SQRT((DBLE(i-vtrack(k,0))**2.0d0) + (DBLE(j-vtrack(k,1))**2.0d0))
                                        !ITS DISTANCE FROM EXISTSING POINT k IS temp
                                        if(temp < 6.0d0) then
                                        !AVERAGING vtrack(k,0),vtrack(k,1) and i,j
                                  		!WRITE (unit=6,fmt="(a, 2f10.2, a , 2i5)") "Averaging ", vtrack(k,0),vtrack(k,1), " and ", i,j
                                        	vtrack(k,0) = (vtrack(k,0)+DBLE(i))/2.0d0
                                    	    vtrack(k,1) = (vtrack(k,1)+DBLE(j))/2.0d0
                                    	!WRITE (unit=6,fmt="(a, 2f10.2)") "We get ", vtrack(k,0),vtrack(k,1)
                                    	end if
                                        if(temp < dist) then
                                        dist = temp
                                        end if
                                end do
                                if(dist > 6.0d0 ) then
                                        vtrack(vlength,0) = DBLE(i)
                                        vtrack(vlength,1) = DBLE(j)
                                        CALL LINEINTVF(velx,vely,i-5,i+5,j-5,j+5,ret)
                             			vtrack(vlength,2) = ret
                                        vlength = vlength + 1
                                      	!WRITE (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                                end if
                        else
                                vtrack(vlength,0) = DBLE(i)
                                vtrack(vlength,1) = DBLE(j)
                                CALL LINEINTVF(velx,vely,i-5,i+5,j-5,j+5,ret)
                                vtrack(vlength,2) = ret
                                vlength = vlength + 1
                                !WRITE (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", uux
                        end if
                end if
        END DO
        END DO
        vort =  vlength
		vortarray = vtrack
        !WRITE (unit=6,fmt="(i5)") vlength
END SUBROUTINE
