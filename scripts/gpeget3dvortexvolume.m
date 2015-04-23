function volume = gpeget3dvortexvolume(dens,gridx,gridy,gridz,th)
    dx = gridx(5)-gridx(4);
    volume = dx*dx*dx*trapz(trapz(trapz(dens<th)))/(range(gridx)*range(gridy)*range(gridz));
end