function h = gpe3disosurface_clean(gridx,gridy,gridz,dens,clean,alphac,th)
dens(clean<th+0.2) = 1.0;
cleanwf = dens;

daspect([1,1,1])
[X,Y,Z] = meshgrid(gridy,gridx,gridz);
p = patch(isosurface(X,Y,Z,cleanwf,th*max(max(max(dens)))));
isonormals(X,Y,Z,cleanwf,p);
set(p,'facecolor','red','edgecolor','none','FaceAlpha',alphac);
camlight;
lighting gouraud;

p = patch(isosurface(X,Y,Z,clean,th*max(max(max(clean)))));
isonormals(X,Y,Z,clean,p);
set(p,'facecolor','blue','edgecolor','none','FaceAlpha',alphac);

axis([-max(gridy) max(gridy) -max(gridx) max(gridx) -max(gridz) max(gridz)]);
view(66, 30);