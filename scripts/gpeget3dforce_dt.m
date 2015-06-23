function [FX,FY,FZ] = gpeget3dforce_dt(dirarg,startno,stride,endno,speed)
for i=startno:stride:endno
    [gridx,gridy,gridz,psi,potential] = gpe3dgetPSI(dirarg,i,speed);
    fprintf('read %d\n',i);
    j = (i+(stride-startno))/stride;
    [fxt,fyt,fzt] = gpeget3dforce(psi,gridx,gridy,gridz,potential);
    FX(j) = fxt;
    FY(j) = fyt;
    FZ(j) = fzt;
end
end