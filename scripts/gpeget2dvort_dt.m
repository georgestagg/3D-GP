function [total,pos,neg] = gpeget2dvort_dt(dirarg,startno,endno,speed,nx,ny)
total = [];
pos = [];
neg = [];
for i=startno:endno
    [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    [xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential);
    total(i+1) = length(pol);
    neg(i+1) = sum(pol);
    pos(i+1) = length(pol) - sum(pol);
end
imagesc(dens)

end