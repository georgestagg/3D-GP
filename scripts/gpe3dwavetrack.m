function [maxn,maxp] = gpe3dwavetrack(dirarg,startno,endno,speed)
     dirarg = regexprep(dirarg, '/$', '');
     pngfolder = strcat(dirarg, '/png');
     for i=startno:endno
        wr = 0;
        [gridx,gridy,gridz,dens,phase,potential] = gpe3dgetWF(dirarg,i,speed);
        for z=1:length(gridz)
            [xlocs,ylocs,pol] = gpeget2dvort_homg(dens(:,:,z),phase(:,:,z),gridx,gridy);
            wr(z) = sqrt(xlocs.^2+ylocs.^2);
        end
        
        hold on
        maxn(i+1) = gridz(find(wr(1:floor(length(wr)/2)) == (max(wr(1:floor(length(wr)/2)))),1));
        maxp(i+1) = gridz(ceil(length(wr)/2)+find(wr(ceil(length(wr)/2):length(wr)) == max(wr(ceil(length(wr)/2):length(wr))),1));
        fprintf('read %d\n',i);
        j = i-startno;
        h=figure('visible','off');
        plot(gridz,wr)
        axis([-max(gridz) max(gridz) 0 1.2]);
        
        filename = strcat(pngfolder, '/w%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        close(h);
     end
end
