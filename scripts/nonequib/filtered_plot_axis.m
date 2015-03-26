%function filtered_plot_axis(i)

clear all;
format long;

Nslice=300;

FileLocation='dens/';
FilePrefix='dens';

M=moviein(Nslice);

for i=1:Nslice
    %open file
    FileName=strcat(FileLocation,FilePrefix,sprintf('%04i',i),'000.dat');
    data=load(FileName);
    
    data_sort=sortrows(data,[1,2,3]);
    
    rangex=floor(min(data_sort(:,1))):1.0:ceil(max(data_sort(:,1)));
    rangey=floor(min(data_sort(:,2))):1.0:ceil(max(data_sort(:,2)));
    rangez=floor(min(data_sort(:,3))):1.0:ceil(max(data_sort(:,3)));
    
    Nx=max(size(rangex));
    Ny=max(size(rangey));
    Nz=max(size(rangez));
    
    n=1;
    for j=1:Nx
        for k=1:Ny
            for l=1:Ny
                u1(j,k,l)=data_sort(n,4)+sqrt(-1)*data_sort(n,5);
                %u2(j,k,l)=data_sort(n,6)+sqrt(-1)*data_sort(n,7);
                n=n+1;
            end
        end
    end
    
    v1=fftn(u1);
    v1=v1/(sqrt(Nx*Ny*Nz));
    v1=fftshift(v1);
    
    %v2=fftn(u2);
    %v2=v2/(sqrt(Nx*Ny*Nz));
    %v2=fftshift(v2);
    
    kx=(-(Nx-1)/2:(Nx-1)/2);
    ky=kx;
    kz=kx;
    
    %cutoff
    kc=8;
    
    for j=1:Nx
        kx2=kx(j)^2;
        for k=1:Ny
            ky2=ky(k)^2;
            for l=1:Ny
                kz2=kz(l)^2;
                k2=kx2+ky2+kz2;
                v1(j,k,l)=v1(j,k,l)*max(1-k2/kc^2,0);
                %v2(j,k,l)=v2(j,k,l)*max(1-k2/kc^2,0);
            end
        end  
    end
    
    v1=ifftshift(v1);
    v1=v1*(sqrt(Nx*Ny*Nz));
    u1=ifftn(v1);
    
    %v2=ifftshift(v2);
    %v2=v2*(sqrt(Nx*Ny*Nz));
    %u2=ifftn(v2);
    
  
    u1=abs(u1).^2;
    %u2=abs(u2).^2;
    
    spatial_av1=mean(mean(mean(u1)));
    value1=0.05*spatial_av1;
    
    %spatial_av2=mean(mean(mean(u2)));
    %value2=0.05*spatial_av2;
    
    p1 = patch(isosurface(rangex,rangey,rangez,u1,value1));
    isonormals(rangex,rangey,rangez,u1,p1);
    axis([-32 32 -32 32 -32 32]);
    set(p1,'facecolor','red','edgecolor','none');
    set(gca,'PlotBoxAspectRatio',[1 1 1],'fontsize',30);
    
    %xlabel('x(\zeta)')
    %ylabel('y(\zeta)')
    %zlabel('z(\zeta)')
    %set(gca,'XTickLabel',[])
    %set(gca,'YTickLabel',[])
    %set(gca,'ZTickLabel',[])
    
    
    %hold on;
    %p2 = patch(isosurface(rangex,rangey,rangez,u2,value2));
    %isonormals(rangex,rangey,rangez,u2,p2);
    %axis([-32 32 -32 32 -32 32]);
    %set(p2,'facecolor','blue','edgecolor','none');
    
    camlight;
    lighting gouraud;
    view(-40,26)
    hold off;
    
    box on;

    M(i)=getframe(gcf);
    clf;
    
end

movie(M);
movie2avi(M, FileSave);



