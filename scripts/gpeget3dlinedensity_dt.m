function volume = gpeget3dlinedensity_dt(dirarg,startno,stride,endno,speed,th)
for i=startno:stride:endno
    [gridx,gridy,gridz,dens,phase,~] = gpe3dgetWF(dirarg,i,speed);
    fprintf('read %d\n',i);
    j = (i+(stride-startno))/stride;
    [ndens, ~] = gpe3disosurface_cutoff(gridx,gridy,gridz,dens,phase,0.8,'red',th,8);
    vol_t = gpeget3dlinedensity(ndens,gridx,gridy,gridz,th);
    volume(j) = vol_t;
end
end