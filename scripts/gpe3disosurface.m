function h = gpe3disosurface(gridx,gridy,gridz,dens,alphac,col,th)
daspect([1,1,1])
[X,Y,Z] = meshgrid(gridy,gridx,gridz);
p = patch(isosurface(X,Y,Z,dens,th*max(max(max(dens)))));
isonormals(X,Y,Z,dens,p);
set(p,'facecolor',col,'edgecolor','none','FaceAlpha',alphac);
camlight;
lighting gouraud;

axis([-max(gridy) max(gridy) -max(gridx) max(gridx) -max(gridz) max(gridz)]);
view(66, 30);