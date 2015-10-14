function dp = gpeget2dvort_disp_image_dt(dirarg,startno,stride,endno,speed)
[~,~,~,dens,~,~] = gpe3dgetWF(dirarg,startno,speed);
dp = trapz(dens,3)/193;
for i=startno+1:stride:endno
    [~,~,~,dens,~,~] = gpe3dgetWF(dirarg,i,speed);
    fprintf('read %d\n',i);    
    dp = (dp + (trapz(dens,3)/193));
end
dp = dp./((endno-startno)/stride);
fclose('all');
end