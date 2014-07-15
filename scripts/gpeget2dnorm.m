function norm = gpeget2dnorm(gridx,gridy,dens)
dx = gridx(2)-gridx(1);
norm = dx*dx*simpson(simpson(dens));
end

