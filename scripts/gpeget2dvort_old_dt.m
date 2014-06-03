function [total,pos,neg] = gpeget2dvort_old_dt(dirarg,startno,endno,speed,nx,ny)
total = [];
pos = [];
neg = [];
for i=startno:endno
    [xlocs,ylocs,pol] = gpeget2dvort_old(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    total(i+1) = length(pol);
    neg(i+1) = -sum(pol-1)/2;
    pos(i+1) = length(pol) - neg(i+1);
end

end