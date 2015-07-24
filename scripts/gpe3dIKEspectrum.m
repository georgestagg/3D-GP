function [kk,Fk] = gpe3dIKEspectrum(gridx,gridy,gridz,psi)

dims = size(psi);
Nx = dims(1);
Ny = dims(2);
Nz = dims(3);
dx=gridx(2)-gridx(1);
dy=gridy(2)-gridy(1);
dz=gridz(2)-gridz(1);
Fk = 0;
kk = 0;
C = {};
[~, v1] = kineticenergydecomposition2d(fftn(psi),dx);

v1=ifftshift(fftn(v1)/(sqrt(Nx*Ny*Nz)));
kx=(-(Nx-1)/2:(Nx-1)/2);
ky=(-(Ny-1)/2:(Ny-1)/2);
kz=(-(Nz-1)/2:(Nz-1)/2);

for m=1:Nx
    for n=1:Ny
        for o=1:Nz
            k=round(sqrt(kx(m)^2+ky(n)^2+kz(o)^2)+0.5);
            if(length(C)<k)
                C{k} = {};
            end
            if(k>2)
                C{k}{end+1} = abs(v1(m,n,o)).^2;
            end
        end
    end
end

for m=1:length(C)
    kk = [kk,m];
    Fk = [Fk,mean(cell2mat(C{m}))];
end