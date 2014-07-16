function [total,pos,neg] = gpeget2dvort_dt(dirarg,startno,stride,endno,speed,nx,ny)
total = [];
pos = [];
neg = [];
for i=startno:stride:endno
    [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    [xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential);
    j = (i+(stride-startno))/stride;
    total(j) = length(pol);
    neg(j) = -sum(pol-1)/2;
    pos(j) = length(pol) - neg(j);
end
%imagesc(dens)
fclose('all');
end