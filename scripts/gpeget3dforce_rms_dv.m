function [FX,FY,FZ] = gpeget3dforce_rms_dv(dirarg,startno,stride,endno)
    datalocation = strcat(dirarg, '%01d');
    for j=1:6
        fname = sprintf(datalocation,j);
        sp =  j*100;
        [ffx,ffy,ffz] = gpeget3dforce_dt(fname,startno,stride,endno,sp);
        FX(j) = mean(ffx)
        FY(j) = mean(ffy)
        FZ(j) = mean(ffz)
    end
end
