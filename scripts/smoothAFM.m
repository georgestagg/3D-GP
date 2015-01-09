x=reshape(x,1024,1024);
y=reshape(y,1024,1024);
z=reshape(z,1024,1024);
h = fspecial('gaussian', size(z), 6);
z2 = imfilter(z, h);
surf(x,y,z2,'EdgeColor', 'none')
x=squeeze(reshape(x,1048576,1));
y=squeeze(reshape(y,1048576,1));
z2=squeeze(reshape(z2,1048576,1));

formatSpec = '%8.6f    %8.6f    %8.6f\n';
fileID = fopen('/data/a8034837/afm-surface/profile1-1um-big-smooth.dat','w');
for i=1:1048576
    fprintf(fileID,formatSpec,x(i),y(i),z2(i));
end
fclose(fileID);