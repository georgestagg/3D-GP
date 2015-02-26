function h = gpe3disosurface(gridx,gridy,gridz,dens,alphac,col)
daspect([1,1,1])
[X,Y,Z] = meshgrid(gridx,gridy,gridz);
p = patch(isosurface(X,Y,Z,dens,0.2*max(max(max(dens)))));
isonormals(X,Y,Z,dens,p);
set(p,'facecolor',col,'edgecolor','none','FaceAlpha',alphac);
camlight;
lighting gouraud;

axis([-max(gridx) max(gridx) -max(gridy) max(gridy) -max(gridz) max(gridz)]);
view(66, 30);