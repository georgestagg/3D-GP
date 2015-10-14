function [gridx,gridy,gridz,dens,phase,potential] = gpe3dgetWF(dirarg,startno,speed)

dirarg = regexprep(dirarg, '/$', '');
if(startno < 10000)
    datalocation = strcat(dirarg, '/%04d.dumpwf.%04d');
    fname = sprintf(datalocation,speed,startno);
elseif (startno < 20000)
    datalocation = strcat(dirarg, '/%04d.dumpwg.%04d');
    fname = sprintf(datalocation,speed,startno-10000);
end

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

dens = psi.*conj(psi);
phase = atan2(imag,real);
phase = permute(phase,[2,1,3]);
potential = permute(potential,[2,1,3]);

fclose('all');
end
