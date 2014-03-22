function [FX,FY] = gpeprint2dforce(dirarg,startno,speed,nx,ny)
dirarg = regexprep(dirarg, '/$', ''); 
datalocation = strcat(dirarg, '/%2dplotwf.%04d');


fname = sprintf(datalocation,speed,startno);
dens = fopen(fname);
A = fscanf(dens, '%g %g %g %g %g\n', [5 inf]);
fclose(dens);
gridy = A(2,1:ny+1);
gridx = A(1,1:ny+1:((ny+1)*(nx+1)));

psi = A(3,:) + 1i*A(4,:);
psi = reshape(psi,ny+1,nx+1);

status ='Loaded frame: %04d...';
status = sprintf(status,startno);
disp(status)
dx = gridx(5)-gridx(4);
dy = gridy(5)-gridy(4);
[FX,FY] = gradient(psi,dx);
FX = conj(psi).*(1i.*FX);
FY = conj(psi).*(1i.*FY);
clear "psi"
FX = trapz(trapz(FX,1),2)*dx*dy;
FY = trapz(trapz(FY,1),2)*dx*dy;

fclose('all');
end