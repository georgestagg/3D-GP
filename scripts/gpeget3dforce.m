function [FX, FY, FZ] = gpeget3dforce(psi,gridx,gridy,gridz,potential)
    dx = gridx(5)-gridx(4);
    dy = gridy(5)-gridy(4);
    dz = gridz(5)-gridz(4);
    dims = size(real(psi));

    [DPsi(1,:,:,:),DPsi(2,:,:,:),DPsi(3,:,:,:)] = gradient(psi,dx,dy,dz);
    [D2Psi(1,1,:,:,:),D2Psi(2,1,:,:,:),D2Psi(3,1,:,:,:)] = gradient(squeeze(DPsi(1,:,:,:)),dx,dy,dz);
    [D2Psi(1,2,:,:,:),D2Psi(2,2,:,:,:),D2Psi(3,2,:,:,:)] = gradient(squeeze(DPsi(2,:,:,:)),dx,dy,dz);
    [D2Psi(1,3,:,:,:),D2Psi(2,3,:,:,:),D2Psi(3,3,:,:,:)] = gradient(squeeze(DPsi(3,:,:,:)),dx,dy,dz);
    
    
    for j=1:3,
    for k=1:3,
        T(j,k,:,:,:) = real((1.0/4.0)*(conj(squeeze(DPsi(j,:,:,:))).*squeeze(DPsi(k,:,:,:))  ...
                                 -conj(psi).*squeeze(D2Psi(j,k,:,:,:))         ...
                                 +squeeze(DPsi(j,:,:,:)).*conj(squeeze(DPsi(k,:,:,:))) ...
                                 -psi.*conj(squeeze(D2Psi(j,k,:,:,:))))        ...
                    +((j==k)/2)*(psi.*conj(psi).*psi.*conj(psi)));
    end
    end

    %imagesc(gridx,gridy,real(squeeze(T(1,2,:,:,30))))


    [DT(1,1,1,:,:,:),~,~] = gradient(squeeze(T(1,1,:,:,:)),dx,dy,dz);
    [DT(1,1,2,:,:,:),~,~] = gradient(squeeze(T(1,2,:,:,:)),dx,dy,dz);
    [DT(1,1,3,:,:,:),~,~] = gradient(squeeze(T(1,3,:,:,:)),dx,dy,dz);
    [~,DT(2,2,1,:,:,:),~] = gradient(squeeze(T(2,1,:,:,:)),dx,dy,dz);
    [~,DT(2,2,2,:,:,:),~] = gradient(squeeze(T(2,2,:,:,:)),dx,dy,dz);
    [~,DT(2,2,3,:,:,:),~] = gradient(squeeze(T(2,3,:,:,:)),dx,dy,dz);
    [~,~,DT(3,3,1,:,:,:)] = gradient(squeeze(T(3,1,:,:,:)),dx,dy,dz);
    [~,~,DT(3,3,2,:,:,:)] = gradient(squeeze(T(3,2,:,:,:)),dx,dy,dz);
    [~,~,DT(3,3,3,:,:,:)] = gradient(squeeze(T(3,3,:,:,:)),dx,dy,dz);

    for k=1:3,
        DsumT(k,:,:,:) = squeeze(DT(1,1,k,:,:,:)+DT(2,2,k,:,:,:)+DT(3,3,k,:,:,:));
    end

    %imagesc(gridx,gridy,squeeze(DsumT(2,:,:,30)));
    [DPot(1,:,:,:),DPot(2,:,:,:),DPot(3,:,:,:)] = gradient(potential,dx,dy,dz);
    
    for k=1:3,
        rhoDv(k,:,:,:) = psi.*conj(psi).*squeeze(DPot(k,:,:,:));
    end
    
    FX = dx*dy*dz*simpson(simpson(simpson(squeeze(DsumT(1,:,:,:)).*(potential>1))))+dx*dy*dz*simpson(simpson(simpson(squeeze(rhoDv(1,:,:,:)))));
    FY = dx*dy*dz*simpson(simpson(simpson(squeeze(DsumT(2,:,:,:)).*(potential>1))))+dx*dy*dz*simpson(simpson(simpson(squeeze(rhoDv(2,:,:,:)))));
    FZ = dx*dy*dz*simpson(simpson(simpson(squeeze(DsumT(3,:,:,:)).*(potential>1))))+dx*dy*dz*simpson(simpson(simpson(squeeze(rhoDv(3,:,:,:)))));
    
end