function dp = gpeget2dvort_displacement_dt(dirarg,startno,stride,endno,speed)
dp = [];
for i=startno:stride:endno
    [gridx,gridy,gridz,dens,phase,~] = gpe3dgetWF(dirarg,i,speed);
    [ndens,nphase] = gpe3disosurface_cutoff(gridx,gridy,gridz,dens,phase,0.8,'red',0.04,10);
    close(gcf);
    fprintf('read %d\n',i);
    disp = mean(gpeget2dvort_displacement(ndens,nphase,gridx,gridy,gridz));
    j = (i+(stride-startno))/stride;
    dp(j) = disp
end
fclose('all');
end