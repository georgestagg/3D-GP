function phase = gpecalcphase(psi)

rpsi = real(psi);
ipsi = imag(psi);

phase = atan2(ipsi./(abs(rpsi)+1e-16));

if(rpsi < 0)
if(ipsi > 0)
    phase = pi - phase;
else
    phase = -pi - phase;
end
end


end