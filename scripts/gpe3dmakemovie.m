function gpe3dmakemovie(dirarg,startno,stride,endno,speed)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,gridz,dens,phase,potential] = gpe3dgetWF(dirarg,i,speed);
        [~,~,~,dens2,~,~] = gpe3dgetWF(dirarg,i+1,speed);
        [~,~,~,dens3,~,~] = gpe3dgetWF(dirarg,i+2,speed);
        [~,~,~,dens4,~,~] = gpe3dgetWF(dirarg,i+3,speed);
        [~,~,~,dens5,~,~] = gpe3dgetWF(dirarg,i+4,speed);
        dens = dens+dens2+dens3+dens4+dens5./5;
        fprintf('read %d\n',i);
        j = i/stride -startno/stride;
        h=figure('visible','off','renderer','painters');
        %gpe3disosurface_cutoff(gridx,gridy,gridz,dens,phase,0.8,'red',0.04,10);
        imagesc(gridx(43:150),gridy(43:150),squeeze(trapz(dens(43:150,43:150,:),3)));
        axis image;
        axis xy;
        filename = strcat(pngfolder, '/p%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
