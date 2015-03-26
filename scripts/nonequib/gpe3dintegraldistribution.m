function [Fx,Fk] = gpe3dintegraldistribution(dens,phase)
Fk=0;
Fx=0;
dims = size(dens);
Nx = dims(1);
Ny = dims(2);
Nz = dims(3);

psi = sqrt(dens).*exp(1i*phase);
v1=fftn(psi);
v1=v1/(sqrt(Nx*Ny*Nz));
v1=ifftshift(v1);
    
kx=(-(Nx-1)/2:(Nx-1)/2);
ky=(-(Ny-1)/2:(Ny-1)/2);
kz=(-(Nz-1)/2:(Nz-1)/2);


for kwv=0:60
    disp(num2str(kwv));
    w1 = v1;
    for j=1:Nx
        for k=1:Ny
            for l=1:Nz
                k2=kx(j)^2+ky(k)^2+kz(l)^2;
                if(k2 > kwv^2)
                    w1(j,k,l)=0;
                end
            end
        end  
    end
    w1=abs(w1).^2;
    Fk(kwv+1)=sum(sum(sum(w1)));
    Fx(kwv+1)=kwv;
end
% figure
% imagesc(squeeze(u1(64,:,:))')
