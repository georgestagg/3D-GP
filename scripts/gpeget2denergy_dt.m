function [totalE,kinE,potE,interE] = gpeget2denergy_dt(dirarg,startno,stride,endno,speed,nx,ny,g)
totalE = [];
kinE = [];
potE = [];
for i=startno:stride:endno
    [gridx,gridy,psi] = gpeget2dPSI(dirarg,i,speed,nx,ny);
    [~,~,~,~,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    [tE,kE,pE,lE] = gpeget2denergy(gridx,gridy,psi,potential,g);
    j = (i+(stride-startno))/stride;
    totalE(j) = tE;
    kinE(j) = kE;
    potE(j) = pE;
    interE(j) = lE;
end
fclose('all');
end