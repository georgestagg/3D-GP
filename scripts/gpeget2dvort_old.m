function [xlocs,ylocs,pol] = gpeget2dvort_old(dirarg,startno,speed,nx,ny)
dirarg = regexprep(dirarg, '/$', ''); 

datalocation = strcat(dirarg, '/%1dplotvort.%04d');
fname = sprintf(datalocation,speed,startno);
densn = fopen(fname);
A = fscanf(densn, '%g %g %g\n', [3 inf]);
A = unique(A','rows');
fclose(densn);

xlocs = A(:,1);
ylocs = A(:,2);
pol =   A(:,3);

fclose('all');