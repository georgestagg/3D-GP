function [totalE,kinE,potE,interE] = gpe3denergy(gridx,gridy,gridz,inpsi,inpot,g)
dx=gridx(2)-gridx(1);
dy=gridy(2)-gridy(1);
dz=gridz(2)-gridz(1);

[FX, FY, FZ] = gradient(inpsi,dx,dy,dz);
dpsi = FX.*conj(FX)+FY.*conj(FY)+FZ.*conj(FZ);
kin=0.5*dpsi;

pot=(inpot).*inpsi.*conj(inpsi);
inter=(g/2).*inpsi.*conj(inpsi).*inpsi.*conj(inpsi);

kinE=kin;
potE=pot;
interE=inter;
totalE=kinE+potE+interE;

end 