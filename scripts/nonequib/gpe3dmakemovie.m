function gpe3dmakemovie(dirarg,startno,stride,endno,speed)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,gridz,dens,phase,potential] = gpe3dgetWF(dirarg,i,speed);
        fprintf('read %d\n',i);
        j = i/stride;
        h=figure('visible','off','renderer','painters');
        gpe3disosurface_cutoff(gridx,gridy,gridz,dens,phase,0.5,'red',0.05,8);
        filename = strcat(pngfolder, '/p%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
