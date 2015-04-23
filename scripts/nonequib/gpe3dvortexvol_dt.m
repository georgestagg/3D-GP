function VX = gpe3dvortexvol_dt(dirarg,startno,stride,endno,speed,th,kc)
for i=startno:stride:endno
    [gridx,gridy,gridz,dens,phase,~] = gpe3dgetWF(dirarg,i,speed);
    [ndens,~] = gpe3d_cutoff(gridx,gridy,gridz,dens,phase,kc);
    fprintf('read %d\n',i);
    j = (i+(stride-startno))/stride;
    VX(j) = gpeget3dvortexvolume(ndens,gridx,gridy,gridz,th)
end
end