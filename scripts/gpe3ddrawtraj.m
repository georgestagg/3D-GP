function gpe3ddrawtraj(gridx,gridy,gridz,phase)

[mgx,mgy,mgz,velx,vely,velz] = gpe3dvelocity(gridx,gridy,gridz,phase);

partx = [1];
party = [1];
partz = [0];

for n=2:10000

   vx = interp3(mgx,mgy,mgz,velx,partx(n-1),party(n-1),partz(n-1));
   vy = interp3(mgx,mgy,mgz,vely,partx(n-1),party(n-1),partz(n-1));
   vz = interp3(mgx,mgy,mgz,velz,partx(n-1),party(n-1),partz(n-1));
   
   partx(n) = partx(n-1) + vx*0.001;
   party(n) = party(n-1) + vy*0.001;
   partz(n) = partz(n-1) + vz*0.001;
end

plot3(partx,party,partz);
axis([-max(gridx) max(gridx) -max(gridy) max(gridy) -max(gridz) max(gridz)]);
end