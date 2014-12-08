function h = gpe3disosurface(gridx,gridy,gridz,dens)
h=figure();
clf;
daspect([1,1,1])
[X,Y,Z] = meshgrid(gridx,gridy,gridz);
p = patch(isosurface(X,Y,Z,dens,0.05*max(max(max(dens)))));
isonormals(X,Y,Z,dens,p);
set(p,'facecolor','red','edgecolor','none');
camlight;
lighting gouraud;
alpha(0.5)
axis([-max(gridx) max(gridx) -max(gridy) max(gridy) -max(gridz) max(gridz)]);
view(66, 30);