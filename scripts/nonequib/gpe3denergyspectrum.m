function [kk,Fk] = gpe3denergyspectrum(dens,phase)

dims = size(dens);
Nx = dims(1);
Ny = dims(2);
Nz = dims(3);
Fk = zeros(1,Nx*Ny*Nz);
kk=[];

psi = sqrt(dens).*exp(1i*phase);
energy = abs(gradient(psi)).^2+real(psi.*conj(psi).*psi.*conj(psi));
v1=fftn(energy);
v1=v1/(sqrt(Nx*Ny*Nz));
v1=ifftshift(v1);
kx=(-(Nx-1)/2:(Nx-1)/2);
ky=(-(Ny-1)/2:(Ny-1)/2);
kz=(-(Nz-1)/2:(Nz-1)/2);

for j=1:Nx
    for k=1:Ny
        for l=1:Nz
            k2=round(kx(j)^2+ky(k)^2+kz(l)^2)+1;
            if(Fk(k2)==0)
                Fk(k2)=(abs(v1(j,k,l)).^2);
            else
                Fk(k2)=(Fk(k2)+abs(v1(j,k,l)).^2)/2.0;
            end
        end
    end
end
end