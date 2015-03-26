function [ndens,nphase] = gpe3disosurface_cutoff(gridx,gridy,gridz,dens,phase,alphac,col,th,kc)
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
%u1=abs(v1).^2;
%imagesc(squeeze(u1(32,:,:))')
for j=1:Nx
    for k=1:Ny
        for l=1:Nz
            k2=kx(j)^2+ky(k)^2+kz(l)^2;
            v1(j,k,l)=v1(j,k,l)*max(1-k2/kc^2,0);
        end
    end  
end
% figure
% u1=abs(v1).^2;
% imagesc(squeeze(u1(32,:,:))')  
v1=ifftshift(v1);
v1=v1*(sqrt(Nx*Ny*Nz));
u1=ifftn(v1);
nphase = atan2(imag(u1),real(u1));
ndens=abs(u1).^2;
spatial_av1=mean(mean(mean(ndens)));
value1=th*spatial_av1;
%figure
daspect([1,1,1]);
p1 = patch(isosurface(gridx,gridy,gridz,ndens,value1));
isonormals(gridx,gridy,gridz,ndens,p1);
set(p1,'facecolor',col,'edgecolor','none','FaceAlpha',alphac);
camlight;
lighting gouraud;
axis([-max(gridy) max(gridy) -max(gridx) max(gridx) -max(gridz) max(gridz)]);
view(66, 30);