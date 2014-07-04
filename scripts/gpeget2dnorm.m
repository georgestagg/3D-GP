function norm = gpeget2dnorm(gridx,gridy,dens)
norm = trapz(gridy,trapz(gridx,dens));
end

