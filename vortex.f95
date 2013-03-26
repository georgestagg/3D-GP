subroutine get_vortices(vort,vortarray)
        use params    
        implicit none
        integer :: i,j,k,vort,vlength
        double precision :: tphase,phase1, phase2, phase3, phase4,temp,dist
        double precision, dimension(0:NX*NY,0:2) :: vtrack,vortarray
        double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase

        vort = 0
        vortarray = 0
        vtrack = 0.0d0
        vlength = 0

        call calc_phase(phase)
       	do i = -NX/2, NX/2-1
		do j = -NY/2, NY/2-1
			call unwrap(phase(i+1,j) - phase(i,j),phase1)
			call unwrap(phase(i+1,j+1) - phase(i+1,j),phase2)
			call unwrap(phase(i,j+1) - phase(i+1,j+1),phase3)
			call unwrap(phase(i,j) - phase(i,j+1),phase4)

			tphase = phase1+phase2+phase3+phase4
			
			if(abs(tphase) > 3.0) then
				write(6,fmt=*) "We have a ", tphase, " at ", i ," , ", j," !"
				write(6,fmt=*) "So...Found one at ", i, ",",j
 	    	   	if(vlength > 0) then
         		   do k = 0, vlength-1
        		    	temp = SQRT((dble(i-vtrack(k,0))**2.0d0) + (dble(j-vtrack(k,1))**2.0d0))
    	    	        !ITS DISTANCE FROM EXISTSING POINT k IS temp
	        	        if(temp*DSPACE < 2.5d0) then
        	    	    	!AVERAGING vtrack(k,0),vtrack(k,1) and i,j
        		        	write (unit=6,fmt="(a, 2f10.2, a , 2i5)") "Averaging ", vtrack(k,0),vtrack(k,1), " and ", i,j
        		        	vtrack(k,0) = (vtrack(k,0)+dble(i))/2.0d0
    	    	        	vtrack(k,1) = (vtrack(k,1)+dble(j))/2.0d0
	        	        	write (unit=6,fmt="(a, 2f10.2)") "We get ", vtrack(k,0),vtrack(k,1)
        	    	    end if
        		        if(temp < dist) then
        	 	           dist = temp
    	     	       end if
	         	    end do

         	   		if(dist*DSPACE > 2.5d0 ) then
         		       vtrack(vlength,0) = dble(i)
        	 	       vtrack(vlength,1) = dble(j)
   		      	       vtrack(vlength,2) = SIGN(1.0d0, tphase)
	         	       vlength = vlength + 1
	         	       write (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", vtrack(vlength,2)
   		      	   end if
   		      	   
	         	else
            	    vtrack(vlength,0) = dble(i)
        	        vtrack(vlength,1) = dble(j)
    	            vtrack(vlength,2) = SIGN(1.0d0, tphase)
	                vlength = vlength + 1
                	write (unit=6,fmt="(i3, a, i3, a, f10.7)") i, " , ", j, " : ", vtrack(vlength,2)
            	end if
            end if
		end do
		end do
		write(6,fmt=*) " done this frame"
		vort =  vlength
		vortarray = vtrack
end subroutine

recursive subroutine unwrap(theta,ret)
	use params
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