function [xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential)
xlocs=[];
ylocs=[];
pol=[];
TFradius = 9.5;
TFpot = (TFradius*TFradius)/2;
dspace=(gridx(2)-gridx(1));
velx(length(dens),length(dens)) = 0;
vely(length(dens),length(dens)) = 0;
presort(length(dens),length(dens)) = 0;
postsort(length(dens),length(dens)) = 0;

for i = 2:length(dens)-1
for j = 2:length(dens)-1
	if (phase(i+1,j)-phase(i-1,j)<-(pi/2.0d0))
		temp1 = phase(i+1,j)-(phase(i-1,j) - pi);
    elseif (phase(i+1,j)-phase(i-1,j)>(pi/2.0d0))
		temp1 = phase(i+1,j)-(phase(i-1,j) + pi);
	else
		temp1 = phase(i+1,j)-phase(i-1,j);
    end
		velx(i,j) = real(temp1)/dspace;
end
end
	
for i = 2:length(dens)-1
for j = 2:length(dens)-1
	if (phase(i,j+1)-phase(i,j-1)<-(pi/2.0d0))
		temp1 = phase(i,j+1)-(phase(i,j-1) - pi);
    elseif (phase(i,j+1)-phase(i,j-1)>(pi/2.0d0))
		temp1 = phase(i,j+1)-(phase(i,j-1) + pi);
	else
		temp1 = phase(i,j+1)-phase(i,j-1);
    end
	vely(i,j) = real(temp1)/dspace;
end
end

for i = 6:length(dens)-6
for j = 6:length(dens)-6
    if(sqrt(gridx(i)*gridx(i)+gridy(j)*gridy(j))<TFradius)
        if((potential(i,j) < TFpot) && (sqrt(velx(i,j)*velx(i,j)+vely(i,j)*vely(i,j))>2.0))
            presort(i,j)=LINEINTVF(velx,vely,i,i+5,j,j+5);
        end
    end
end
end


%densblob=dens>0.002;
  
    
%phase(dens<0.0005)=0;
%imagesc(sqrt(velx.*velx+vely.*vely))
%imagesc(gridx,gridy,dens)
imagesc(dens)
%imagesc(potential)
%imagesc(presort)

function ret = LINEINTVF(fieldx,fieldy,x,ex,y,ey)
	l1=0.0d0;
	l2=0.0d0;
	l3=0.0d0;
	l4=0.0d0;
	for t = y:ey
		l1 = l1 + dspace*fieldy(x,t);
    end
	for t = x:ex
		l2 = l2 + dspace*fieldx(t,y);
    end
	for t = y:ey
		l3 = l3 + dspace*fieldy(ex,t);
    end
	for t = x:ex
		l4 = l4 + dspace*fieldx(t,ey);
    end
	ret = l2+l3-l4-l1;
end


end

