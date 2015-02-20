function params(m,as,N,wx,wy,wz,r_0) %atom mass, scattering length, # atoms, xyz trap freq, r_0!=0 for ring traps
    format LONGG
    k_B = 1.3806488e-23; 
    hbar=1.0546e-34;
    %m=22.98977*1.6605e-27          % Atomic mass of sodium
    %m=86.90918*1.6605e-27          % Atomic mass of rubidium87
    %as=2.8e-9;                     % Scattering length of sodium   
    %as=96*0.529e-10;               % Scattering length of rubidium87
    
    g_3d = 4*pi*hbar*hbar*as/m;  % Calculate g
    ombar = (wx*wy*wz)^(1/3);
    wx_h0 = wx/ombar
    wy_h0 = wy/ombar
    wz_h0 = wz/ombar   
    if(r_0 > 1e-9)
       mu_3d = hbar*ombar*sqrt(2*N*as/(pi*r_0));
    else
       mu_3d = (hbar*ombar/2)*(15*N*as/sqrt(hbar/(m*ombar)))^(2/5);
    end  
    l_x = sqrt(hbar/(m*wx));
    l_y = sqrt(hbar/(m*wy));
    l_z = sqrt(hbar/(m*wz));
    lbar = (l_x*l_y*l_z)^(1/3)
    
    g_3d_h0 = g_3d*N*(lbar^-3)/(hbar*ombar)
    mu_3d_ho = mu_3d/(hbar*ombar)
    
    g_2d = g_3d/(sqrt(2*pi)*l_z);
    mu_2d = mu_3d - hbar*wz/2;
    g_2d_ho = N*g_2d/(hbar*wx*l_x*l_y)
    mu_2d_ho = mu_2d/(hbar*wx)
    gamma = 0.001;
    T = (gamma*pi*hbar^2)/(3*4*m*as^2*k_B*ombar)*1e9
end