function gpe2ddensvort(dirarg,startno,speed,nx,ny)
clf;
[gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,startno,speed,nx,ny);
[xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential);
imagesc(gridx,gridy,dens)
colormap 'gray'
axis image
axis xy
hold on
g = gscatter(xlocs+0.2,ylocs+0.1,pol,['b','r'],['^','o'],5,'off');
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
axis([-20 20 -20 20]);
axis off
curtime = num2str(roundn(((startno*200)*0.0005)/(15*2*pi)*1000,-1));
text(0,-18,strcat('\color{white}',curtime,' ms'),'FontSize',16,'HorizontalAlignment','center')
end