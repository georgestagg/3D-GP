subroutine calc_OBJPOT (obj)				
	use params
	implicit none
	integer :: i,j,obj,slicei,slicej
	double precision :: di,dj,xi1,ity,xtt,ytt,ztt
	double precision, dimension(1:256) :: xdat,ydat
	double precision, dimension(256,256) :: ydatb
	open (unit = 5, file = "../profile1-1um.dat")
	
	xi1 = 0.066d0 !nanometres
	
	do j = 1,256
	do i = 1,256
		read (5,*) xtt,ytt,ztt
			xdat(i) = ((xtt-0.5d0)*35.0d0)/xi1
			ydatb(i,j) = (ztt/xi1)   !level y so that data=0 -> -DSPACE*NY/2
	end do
	end do
	!ydat = ydatb(slicei,slicej)
	ydat = ydatb(120,:)
	ydat= ydat - MINVAL(ydat(3:253)) - (NY/2)*DSPACE
	!write(6,*) ydat
	OBJPOT = 0.0d0
	
	do j = -NY/2,NY/2
		dj = j*DSPACE
		do i = -NX/2,NX/2
			di = i*DSPACE
			if(di>maxval(xdat(1:252))) then
				di = xdat(252)
			end if
			CALL interp(di,xdat,ydat,ity)
			if(dj < ity) then			
				OBJPOT(i,j) = 50.0d0
			end if
		end do
	end do
	close(5)
end subroutine


subroutine interp(xloc,xdat,ydat,rety)
	implicit none
	integer :: i,si
	double precision :: xloc,rety,low,high,diff,alpha
	double precision, dimension(1:255) :: xdat,ydat
	
	si = 1
	do i = 1,255
		if(xdat(i)>xloc) then
			si = i
			exit
		end if
	end do
	
	low = xdat(si-1)
	high = xdat(si)
	diff = high-low
	alpha = (xloc-low)/diff
	low = ydat(si-1)
	high = ydat(si)
	rety = low+((high-low)*alpha)
end subroutine

