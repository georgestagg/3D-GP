function [gridx,gridy,gridz,psi,potential] = gpe3dgetPSI(dirarg,startno,speed)

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
psi = permute(psi,[2,1,3]);
potential = permute(potential,[2,1,3]);

fclose('all');
end
