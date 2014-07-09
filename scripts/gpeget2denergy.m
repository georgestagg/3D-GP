function [totalE,kinE,potE,lastE] = gpeget2denergy(gridx,gridy,inpsi,inpot)
dx=gridx(2)-gridx(1);
[FX, FY] = derivative5(inpsi, 'x', 'y');
FX = FX/dx;
FY = FY/dx;
%[FX, FY] = gradient(inpsi,dx,dx);
dpsi = FX.*conj(FX)+FY.*conj(FY);
kin=0.5*dpsi;

%[FX, FY] = derivative5(sqrt(inpsi.*conj(inpsi)), 'x', 'y');
%FX = FX/dx;
%FY = FY/dx;
%qp=0.5*(FX.*conj(FX)+FY.*conj(FY));
%qpE=dx*dx*simpson(simpson(qp))

pot=(inpot).*inpsi.*conj(inpsi);
last=(22543.29/2).*inpsi.*conj(inpsi).*inpsi.*conj(inpsi);

kinE=dx*dx*simpson(simpson(kin));
potE=dx*dx*simpson(simpson(pot));
lastE=dx*dx*simpson(simpson(last));
totalE=kinE+potE+lastE;

end 