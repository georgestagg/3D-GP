function gpe2dmakemovie(dirarg,startno,stride,endno,speed,nx,ny)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
        %[xlocs,ylocs,pol] = gpeget2dvort_homg(dens,phase,gridx,gridy,potential);
        fprintf('read %d\n',i);
        j = i/stride;
        h=figure('visible','off');
        imagesc(gridx(1:nx+1),gridy(1:ny+1),dens);
        %colormap(hot);
        
        x=linspace(0,1,128);
        r = sqrt(x);
        g = x.^3.0;
        b=sin(x*2*pi);
        b(b<0)=0;
        pm3d7515_bbry=[r;g;b]';
        colormap(pm3d7515_bbry);
        caxis([0,1])
        axis image;
        axis xy;
        hold on;
        %g = gscatter(xlocs,ylocs-gridy(1),pol,['b','r'],['^','o'],5,'off');
        if(length(g)==1)
	 set(g(1), 'MarkerFaceColor', 'r')
         set(g(1),'Marker','o');
	 set(g(1),'MarkerEdgeColor','none');
	end
	if(length(g)==2)
	  set(g(1),'MarkerEdgeColor','none');
	  set(g(1), 'MarkerFaceColor', 'b')
	  set(g(2),'MarkerEdgeColor','none');
	  set(g(2), 'MarkerFaceColor', 'r')
    end
        xlabel('x/\xi', 'FontSize',16);
        ylabel('y/\xi', 'FontSize',16);
        axis([gridx(1) gridx(nx+1) gridy(1) gridy(ny+1)]);
        filename = strcat(pngfolder, '/p%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
