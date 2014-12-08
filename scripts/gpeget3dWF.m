function [gridx,gridy,gridz,dens,phase,potential] = gpeget3dWF(dirarg,startno,speed,nx,ny,nz)
gridx=0;
gridy=0;
gridz=0;
FX = 0;
FY = 0;
FZ = 0;
dirarg = regexprep(dirarg, '/$', '');

datalocation = strcat(dirarg, '/%04d.dumpwf.%04d');
fname = sprintf(datalocation,speed,startno);
densn = fopen(fname);
A = fscanf(densn, '%g %g %g %g %g %g\n', [6 inf]);
fclose(densn);
gridz = A(3,1:nz+1);
gridy = A(2,1:nz+1:((ny+1)*(nz+1)));
gridx = A(1,1:((ny+1)*(nz+1)):((nx+1)*(ny+1)*(nz+1)));

psi = A(4,:) + 1i*A(5,:);
psi = permute(reshape(psi,(nz+1),(ny+1),(nx+1)),[2,3,1]);
dens = psi.*conj(psi);
phase = atan2(imag(psi),real(psi));
potential = permute(reshape(A(6,:),(nz+1),(ny+1),(nx+1)),[2,3,1]);

fclose('all');
end
