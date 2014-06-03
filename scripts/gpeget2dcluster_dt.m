function [clus_dt] = gpeget2dcluster_dt(dirarg,startno,endno,speed,nx,ny)

for i=startno:10:endno
    [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    k = i/10;
    [xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential);
    [clusters] = gpeget2dcluster(xlocs,ylocs,pol);
    clus_dt(k) = 0;
    for j=1:length(clusters)
        clus_dt(k) = clus_dt(k) + length(clusters{j});
    end
    clus_dt(k) = clus_dt(k)/length(clusters);
end

end