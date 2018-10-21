%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Dhanush VITTAL SHENOY
%               Last modified
%               21/Oct/2018
%
%        https://github.com/Dhanushenoy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
close all;
clf;

%computational parameters------------------------ 
Lx=0.7; Ly=0.1; tmax=1; nu=0.1; dt=0.0001; Uwall=1;
nx=71; ny=11;
%------------------------------------------------

%Raynolds number
Re=1/nu;

%Iteration parameter-----------------------------
maxIt=1; Beta=1.5; maxError=0.001; 
%------------------------------------------------

%----Initialize----------------------------------
strf=zeros(nx,ny);
vort=zeros(nx,ny);
strCalc=zeros(nx,ny);
u=zeros(nx,ny);
v=zeros(nx,ny);
dx=Lx/(nx-1); dy=Ly/(ny-1);
t=0.0;
counter=0;

%For mesh generation
for i=1:nx;for j=1:ny;
        xi(i)=Lx*(i-1)/(nx);
        yi(j)=Ly*(j-1)/(ny);
end;end
[x y]=meshgrid(xi,yi);

%For animation
p=tmax;
F(p/0.01) = struct('cdata',[],'colormap',[]);

%initialize streamfunction
for i=1:nx; for j=1:ny;
   strf(i,j)=((1/ny)*(j));
end;end;


for t=0:dt:tmax
    
%%%%%%%%%%%%%%%Beginning of SOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for It=1:maxIt
    
 strCalc=strf;

 %-------------------Laplacian iteration-------------------------
    for i=2:nx-1;
        for j=2:ny-1
            
            strf(i,j)=0.25*Beta*(strf(i+1,j)+strf(i-1,j)+strf(i,j+1)...
                +strf(i,j-1)+dx*dx*vort(i,j))+(1.-Beta)*strf(i,j);
        end
    end

%------------------Convergece calcualtion------------------------    
    err=0.0;
    for i=1:nx; for j=1:ny;
            err=err+abs(strCalc(i,j)-strf(i,j));
        end;end    
    if err<=maxError;
        break
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------Vorticity boundary for Channel----------------------------
%Top and bottom wall
vort(:,1) = 2*(strf(:,1) - strf(:,2))/(dx*dx);  %Bottom wall
vort(:,ny) = 2*(strf(:,ny) - strf(:,ny-1))/(dx*dx); %Top wall

%Inlet interpolation
 for j=2:ny-1;
    vort(1,2:ny-1)=(strf(1,1:ny-2) + strf(1,3:ny) - 2*strf(1,2:ny-1))/(dy*dy) + 2*(strf(1, 2:ny-1) - strf(2, 2:ny-1))/(dy*dy);
 end

 %Simple outlet [*Need to change*]
vort(nx,2:ny-1)=0;
%----------------------------------------------------------------


%------------------Vorticity calculation-------------------------
for i=2:nx-1; for j=2:ny-1;
        vortCalc(i,j)=((strf(i,j+1)-strf(i,j-1))*(vort(i+1,j)-vort(i-1,j))-...
            (strf(i+1,j)-strf(i-1,j))*(vort(i,j+1)-vort(i,j-1)))/(4.0*dx*dx)+...
        (1/Re)*(vort(i+1,j)+vort(i-1,j)+vort(i,j+1)+vort(i,j-1)-4.0*vort(i,j))/(dx*dx);
end;end;

%Time march--------------------------------------------------------
vort(2:nx-1,2:ny-1)=vort(2:nx-1,2:ny-1)+dt*vortCalc(2:nx-1,2:ny-1);
t=t+dt;
%------------------------------------------------------------------

%To print------------
counter=t/tmax*100;
display(counter,'%')
Re
%--------------------

if mod(t,0.001) == 0
%For animation
contourf(x,y,rot90(fliplr(vort)),'EdgeColor','none');
axis equal
colormap('jet')
colorbar;
drawnow
end
end



%--------------velocity----------------
%
%   Velocity from stream function
%
%--------------------------------------
for i=2:nx-1;
    for j=2:ny-1;
    u(i,j)=(strf(i,j+1)-strf(i,j-1))/(dy);
    v(i,j)=-(strf(i+1,j)-strf(i-1,j))/(dx);
    
    end
end
magV=sqrt(u.^2+v.^2);   %magnitude for plots
%--------------------------------------

figure('rend','painters','pos',[10 10 1400 200])
contour(x,y,rot90(fliplr(strf)),'EdgeColor','none');
title('StreamFunction, $R_e=10$','interpreter','latex')
