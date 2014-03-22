function [FX,FY] = force_2d_dt(loc,speed,endno,jump,nx,ny,dt)

for i=jump:jump:endno
        [dtx,dty] = gpeprint2dforce(loc,i,speed,nx,ny);
        FX((i/jump),(j/sjump)+1,(k/sjump)+1) = dtx;
        FY((i/jump),(j/sjump)+1,(k/sjump)+1) = dty;
end
FX = real(diff(FX))/(2*jump*dt);
FY = real(diff(FY))/(2*jump*dt);
end