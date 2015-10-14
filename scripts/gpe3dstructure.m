function [rcr,ducr] = gpe3dstructure(gridx,gridy,gridz,velx,vely,velz)
    dx = gridx(2)-gridx(1);
    dy = gridy(2)-gridy(1);
    dz = gridz(2)-gridz(1);
    du = [];
    ducr = [];
    rcr = [];
    r = [];
    order = 3;
    maxrpx = 16;
    step = 1;
    
    disp('Setup complete. Interating over grid...')
    
    for i=1:5:length(gridx)-2;
    disp(['x=',int2str(i)]);
    drawnow('update');
    for j=1:5:length(gridy)-2;
    for k=1:5:length(gridz)-2;
        for a=max(i-maxrpx,1):step:min(i+maxrpx,length(gridx)-2);
        for b=max(j-maxrpx,1):step:min(j+maxrpx,length(gridy)-2);
        %for c=max(k-maxrpx,1):step:min(k+maxrpx,length(gridz)-2);
         c=k;   
            l = [(a-i)*dx,(b-j)*dy,(c-k)*dz];

            u_x  = [velx(i,j,k),vely(i,j,k),velz(i,j,k)];
            u_xl = [velx(a,b,c),vely(a,b,c),velz(a,b,c)];
            lhat = l'/sqrt(l*l');
            du(end+1) = ((u_xl - u_x)*lhat).^order;
            r(end+1)  = sqrt(l*l');
        %end
        end
        end
    end
    end
    end
    disp('Soring by radius of step vectors...')
    for run=unique(r)
        ducr(end+1) = mean(du(r==run));
        rcr(end+1) = run;
    end  

end
