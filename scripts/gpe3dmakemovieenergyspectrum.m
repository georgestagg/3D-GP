function gpe3dmakemovieenergyspectrum(dirarg,startno,stride,endno,speed)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,gridz,psi,potential] = gpe3dgetPSI(dirarg,i,speed);
        fprintf('read %d\n',i);
        j = i/stride -startno/stride;
        h=figure('visible','off','renderer','painters');
        [kk,Fk] = gpe3denergyspectrum(gridx,gridy,gridz,psi,potential);
        xlabel('k', 'FontSize',16);
        ylabel('P', 'FontSize',16);
        loglog(kk(1:334),tsmovavg(Fk(1:334),'simple',10),'LineWidth',2);
        hold on
        axis([4 170 0.01 25]);
        plot(10:100,300*(10:100).^(-1),'k--')
        plot(10:100,700*(10:100).^(-5/3),'k--')
        text(120,700*100.^(-5/3),'{\it k}^{-5/3}','FontSize',16,'HorizontalAlignment','center','Interpreter','tex')
        text(120,300*100.^(-1),'{\it k}^{-1}','FontSize',16,'HorizontalAlignment','center','Interpreter','tex')
        axes('Position',[.15 .15 .35 .35])
        gpe3disosurface_cutoff(gridx,gridy,gridz,abs(psi).^2,atan2(imag(psi),real(psi)),0.8,'red',0.05,8);
        filename = strcat(pngfolder, '/p%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
