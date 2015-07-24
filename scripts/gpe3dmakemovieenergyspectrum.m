function gpe3dmakemovieenergyspectrum(dirarg,startno,stride,endno,speed)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,gridz,psi,potential] = gpe3dgetPSI(dirarg,i,speed);
        fprintf('read %d\n',i);
        j = i/stride -startno/stride;
        [ndens, ~] = gpe3disosurface_cutoff(gridx,gridy,gridz,abs(psi).^2,atan2(imag(psi),real(psi)),0.8,'red',0.05,8);
        h=figure('visible','off','renderer','painters');
        [kk,Fk] = gpe3denergyspectrum(gridx,gridy,gridz,psi,potential);
        xlabel('k', 'FontSize',16);
        ylabel('P', 'FontSize',16);
        loglog(kk,tsmovavg(Fk,'simple',5),'LineWidth',2);
        hold on
        axis([0 181 0.02 11]);
        plot(10:300,200*(10:300).^(-5/3),'k--')
        %text(120,2000*300.^(-5/3),'{\it k}^{-5/3}','FontSize',16,'HorizontalAlignment','center','Interpreter','tex')
        lined = gpeget3dlinedensity(ndens,gridx,gridy,gridz,0.04);
        kl = (2*pi)./(1./sqrt(lined));
        plot([kl,kl],[0,11],'k--')
        axes('Position',[.15 .15 .35 .35]);
        gpe3disosurface_cutoff(gridx,gridy,gridz,abs(psi).^2,atan2(imag(psi),real(psi)),0.8,'red',0.05,8);
        filename = strcat(pngfolder, '/p%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
