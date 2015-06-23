function [kk,Fk] = gpe3denergyspectrum(gridx,gridy,gridz,psi,potential)

dims = size(psi);
Nx = dims(1);
Ny = dims(2);
Nz = dims(3);
Fk = zeros(1,27649);
kk = [];
dx=gridx(2)-gridx(1);
dy=gridy(2)-gridy(1);
dz=gridz(2)-gridz(1);

[FX, FY, FZ] = gradient(psi,dx,dy,dz);
dpsi = FX.*conj(FX)+FY.*conj(FY)+FZ.*conj(FZ);
kin=0.5*dpsi;

pot=(potential).*psi.*conj(psi);
inter=(1/4).*psi.*conj(psi).*psi.*conj(psi);

kinE=kin;
potE=pot;
interE=inter;
totalE=kinE;%+potE+interE;

%energy = abs(gradient(psi)).^2+real(psi.*conj(psi).*psi.*conj(psi));
v1=fftn(totalE);
v1=v1/(sqrt(Nx*Ny*Nz));
v1=ifftshift(v1);
kx=(-(Nx-1)/2:(Nx-1)/2);
ky=(-(Ny-1)/2:(Ny-1)/2);
kz=(-(Nz-1)/2:(Nz-1)/2);

for j=1:Nx
    for k=1:Ny
        for l=1:Nz
            k2=round(sqrt(kx(j)^2+ky(k)^2+kz(l)^2)*2.0)+1;
            if(k2>6)
                kk(k2) = 1;
                if(Fk(k2) > 0) 
                    Fk(k2) = (Fk(k2)+abs(v1(j,k,l)).^2)/2.0;
                else
                    Fk(k2) = abs(v1(j,k,l)).^2;
                end
            end
        end
    end
end

kk = ((1:length(kk))-1)*0.5;
end