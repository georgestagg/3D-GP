function [Ekinsq_c, Ekinsq_i] = kineticenergydecomposition2d(psik,dx)
% psik - condensate wavefunction in space dimensions [xx,yy,zz]
% dx - lattice spacing - assumes dx = dy = dz and dimensions x = y/z
% outputs: Ekinsq_c, Ekinsq_i
% where:
% Ekinsq_c: Compressible kinetic energy density dimensions [mm,mm,mm]
% Ekinsq_i: Incompressible kinetic energy density dimensions [mm,mm,mm]

gridDim = size(psik);
mm = gridDim(1);

%--------------------------------------------------------------------------
% A: Calculating Velocities
%--------------------------------------------------------------------------

dk = 2*pi/(dx*mm);

ki = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];
kj = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];
kk = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];

psi2 = ifftn(psik);
conjpsi2 = conj(psi2);
dens2 = psi2.*conjpsi2;

psiconjk = (fftn(conjpsi2));

kxpsi= zeros(mm,mm,mm);
kxpsiconj = zeros(mm,mm,mm);
kypsi = zeros(mm,mm,mm);
kypsiconj = zeros(mm,mm,mm);
kzpsi = zeros(mm,mm,mm);
kzpsiconj = zeros(mm,mm,mm);

for ii = 1:1:length(ki)
  for jj = 1:1:length(kj)
        for ll = 1:1:length(kk)
          kxpsi(ii,jj,ll)=(1i*ki(ii))*psik(ii,jj,ll);
          kxpsiconj(ii,jj,ll)=(1i*ki(ii))*psiconjk(ii,jj,ll);
          kypsi(ii,jj,ll)=(1i*kj(jj))*psik(ii,jj,ll);
          kypsiconj(ii,jj,ll)=(1i*kj(jj))*psiconjk(ii,jj,ll);
          kzpsi(ii,jj,ll)=(1i*kk(ll))*psik(ii,jj,ll);
          kzpsiconj(ii,jj,ll)=(1i*kk(ll))*psiconjk(ii,jj,ll);  
        end
  end
end

velx = real(-0.5*1i*(conjpsi2(:,:,:).*ifftn(kxpsi)-psi2(:,:,:).*ifftn(kxpsiconj))./dens2(:,:,:));
vely = real(-0.5*1i*(conjpsi2(:,:,:).*ifftn(kypsi)-psi2(:,:,:).*ifftn(kypsiconj))./dens2(:,:,:));
velz = real(-0.5*1i*(conjpsi2(:,:,:).*ifftn(kzpsi)-psi2(:,:,:).*ifftn(kzpsiconj))./dens2(:,:,:));

clear kxpsi; clear kypsi; clear kzpsi;
clear kxpsiconj; clear kypsiconj; clear kzpsiconj;
clear psik; clear conjpsik; 
clear conjpsi2; clear psi2; 

%--------------------------------------------------------------------------
% B: Calculating the incompressible and compressible components
%--------------------------------------------------------------------------

omegax = sqrt(dens2(:,:,:)).*(velx(:,:,:));
omegay = sqrt(dens2(:,:,:)).*(vely(:,:,:));
omegaz = sqrt(dens2(:,:,:)).*(velz(:,:,:));

clear dens2; 

omegax_kx = (fftn(omegax));
omegay_ky = (fftn(omegay));
omegaz_kz = (fftn(omegaz));

komegac_kx = zeros(mm,mm,mm);
komegac_ky = zeros(mm,mm,mm);
komegac_kz = zeros(mm,mm,mm);

komegai_kx = zeros(mm,mm,mm);
komegai_ky = zeros(mm,mm,mm);
komegai_kz = zeros(mm,mm,mm);

absk = zeros(mm,mm);

for ii = 1:1:length(ki)
  for jj = 1:1:length(kj)
       for ll = 1:1:length(kk)
          absk(ii,jj,ll) = ki(ii)*ki(ii)+kj(jj)*kj(jj)+kk(ll)*kk(ll);
          komegac_ky(ii,jj,ll) = (kj(jj)*ki(ii)*omegax_kx(ii,jj,ll)+kj(jj)*kj(jj)*omegay_ky(ii,jj,ll)+kj(jj)*kk(ll)*omegaz_kz(ii,jj,ll))/(absk(ii,jj,ll));
          komegac_kx(ii,jj,ll) = (ki(ii)*ki(ii)*omegax_kx(ii,jj,ll)+ki(ii)*kj(jj)*omegay_ky(ii,jj,ll)+ki(ii)*kk(ll)*omegaz_kz(ii,jj,ll))/(absk(ii,jj,ll));
          komegac_kz(ii,jj,ll) = (kk(ll)*ki(ii)*omegax_kx(ii,jj,ll)+kk(ll)*kj(jj)*omegay_ky(ii,jj,ll)+kk(ll)*kk(ll)*omegaz_kz(ii,jj,ll))/(absk(ii,jj,ll));
          komegai_ky(ii,jj,ll) = omegay_ky(ii,jj,ll) - komegac_ky(ii,jj,ll);
          komegai_kz(ii,jj,ll) = omegaz_kz(ii,jj,ll) - komegac_kz(ii,jj,ll);
          komegai_kx(ii,jj,ll) = omegax_kx(ii,jj,ll) - komegac_kx(ii,jj,ll);
       end
   end
end

komegac_kz(find(isnan(komegac_kz))) = 0;
komegac_ky(find(isnan(komegac_ky))) = 0;
komegac_kx(find(isnan(komegac_kx))) = 0;

komegai_kz(find(isnan(komegai_kz))) = 0;
komegai_ky(find(isnan(komegai_ky))) = 0;
komegai_kx(find(isnan(komegai_kx))) = 0;

omegac_x = real(ifftn(komegac_kx));
omegac_y = real(ifftn(komegac_ky));
omegac_z = real(ifftn(komegac_kz));

omegai_x = real(ifftn(komegai_kx));
omegai_y = real(ifftn(komegai_ky));
omegai_z = real(ifftn(komegai_kz));

Ekinsq_c = 0.5*((omegac_x.^2+omegac_y.^2+omegac_z.^2));
Ekinsq_i = 0.5*((omegai_x.^2+omegai_y.^2+omegai_z.^2));
end