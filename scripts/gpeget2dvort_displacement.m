function disp = gpeget2dvort_displacement(dens,ophase,gridx,gridy,gridz)
for i=1:length(gridz)
[xl,yl,~] = gpeget2dvort_homg(dens(:,:,i),ophase(:,:,i),gridx,gridy);
disp(i) = sqrt(min(abs(xl))^2 + min(abs(yl))^2);
end
disp = rms(disp - mean(disp));
end

