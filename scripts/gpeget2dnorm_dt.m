function norm = gpeget2dnorm_dt(dirarg,startno,stride,endno,speed,nx,ny)
norm = [];
for i=startno:stride:endno
    [gridx,gridy,dens,~,~] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    normt = gpeget2dnorm(gridx,gridy,dens);
    j = i/stride;
    norm(j) = normt;
end
fclose('all');
end