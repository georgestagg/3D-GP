function [gridx,gridy,gridz,dens,phase,potential] = gpe3dgetWF(dirarg,startno,speed)

dirarg = regexprep(dirarg, '/$', '');
datalocation = strcat(dirarg, '/%04d.dumpwf.%04d');
fname = sprintf(datalocation,speed,startno);

%i=netcdf.open(fname);
gridx = ncread(fname,'x');
gridy = ncread(fname,'y');
gridz = ncread(fname,'z');
real = ncread(fname,'real');
imag = ncread(fname,'imag');
potential = ncread(fname,'pot');
psi = real + 1i.*imag;

%psi = A(4,:) + 1i*A(5,:);
%psi = permute(reshape(psi,(nz+1),(ny+1),(nx+1)),[2,3,1]);

dens = psi.*conj(psi);
phase = atan2(imag,real);
%potential = permute(reshape(A(6,:),(nz+1),(ny+1),(nx+1)),[2,3,1]);

fclose('all');
end
