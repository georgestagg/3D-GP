function [FX, FY,FZ] = gpeget2dforce(dens,phase,gridx,gridy,gridz,potential)
    psi = sqrt(dens) + exp(1i*phase);
    dx = gridx(5)-gridx(4);
    dims = size(dens);

    [DPsi(1,:,:,:),DPsi(2,:,:,:),DPsi(3,:,:,:)] = gradient(psi,dx);

    [D2Psi(1,1,:,:,:),D2Psi(2,1,:,:,:),D2Psi(3,1,:,:,:)] = gradient(DPsi(1,:,:,:),dx);
    [D2Psi(1,2,:,:,:),D2Psi(2,2,:,:,:),D2Psi(3,2,:,:,:)] = gradient(DPsi(2,:,:,:),dx);
    [D2Psi(1,3,:,:,:),D2Psi(2,3,:,:,:),D2Psi(3,3,:,:,:)] = gradient(DPsi(3,:,:,:),dx);
    
    for j=1:3,
    for k=1:3,
        T(j,k,:,:,:) = real((1.0/4.0)*(conj(squeeze(DPsi(j,:,:,:))).*squeeze(DPsi(k,:,:,:))  ...
                                 -conj(psi(:,:,:)).*squeeze(D2Psi(j,k,:,:,:))         ...
                                 +squeeze(DPsi(j,:,:,:)).*conj(squeeze(DPsi(k,:,:,:)))...
                                 -psi(:,:,:).*conj(squeeze(D2Psi(j,k,:,:,:))))        ...
                    +((j==k)/2)*(abs(psi(:,:,:)).^4));
    end
    end


    [DT(1,1,1,:,:,:),~,~] = gradient(squeeze(T(1,1,:,:,:)),dx);
    [DT(1,1,2,:,:,:),~,~] = gradient(squeeze(T(1,2,:,:,:)),dx);
    [DT(1,1,3,:,:,:),~,~] = gradient(squeeze(T(1,3,:,:,:)),dx);
    [~,DT(2,2,1,:,:,:),~] = gradient(squeeze(T(2,1,:,:,:)),dx);
    [~,DT(2,2,2,:,:,:),~] = gradient(squeeze(T(2,2,:,:,:)),dx);
    [~,DT(2,2,3,:,:,:),~] = gradient(squeeze(T(2,3,:,:,:)),dx);
    [~,~,DT(3,3,1,:,:,:)] = gradient(squeeze(T(3,1,:,:,:)),dx);
    [~,~,DT(3,3,2,:,:,:)] = gradient(squeeze(T(3,2,:,:,:)),dx);
    [~,~,DT(3,3,3,:,:,:)] = gradient(squeeze(T(3,3,:,:,:)),dx);
    for k=1:3,
        DsumT(k,:,:,:) = squeeze(DT(1,1,k,:,:,:)+DT(2,2,k,:,:,:)+DT(3,3,k,:,:,:));
    end

    imagesc(gridx,gridy,squeeze(DsumT(1,:,:,128)));
    
    [DPot(1,:,:,:),DPot(2,:,:,:),DPot(3,:,:,:)] = gradient(potential,dx);
    
    for k=1:3,
        rhoDv(k,:,:,:) = dens.*squeeze(DPot(k,:,:,:));
    end 
    %imagesc(gridx,gridy,squeeze(DsumT(1,:,:,:))+squeeze(rhoDv(1,:,:,:)));
    
    FX = -dx*dx*dx*trapz(trapz(trapz(squeeze(DsumT(1,:,:,:)))))-dx*dx*dx*trapz(trapz(trapz(squeeze(rhoDv(1,:,:,:)))));
    FY = -dx*dx*dx*trapz(trapz(trapz(squeeze(DsumT(2,:,:,:)))))-dx*dx*dx*trapz(trapz(trapz(squeeze(rhoDv(2,:,:,:)))));
    FZ = -dx*dx*dx*trapz(trapz(trapz(squeeze(DsumT(3,:,:,:)))))-dx*dx*dx*trapz(trapz(trapz(squeeze(rhoDv(3,:,:,:)))));
end