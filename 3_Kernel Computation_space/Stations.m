%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             CREATE STATIONS (FILE) FOR SPECFEM2D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Useful when you want to plant many receivers or in an specific pattern

close all
clear all
% {
%% -----------------Kernel computation in space---------------------------------
%------------------------------------------
% Parameters - Kernel computation in space
nx=865%216 %410;      %number GLL points in the model (grids!!)
L=16800; %9000;      %length of the medium in meters. We assumed square ->A=LxL   
nlines=70%100;   %number of RCV Lines (#1 at the top, last one at the bottom)
nrec=70%100;     %number of rcv for each line (#1 at the LHS)
pml=0%4;  %PML boundary conditions in terms of number of spectral elements nx. If no borders write 0

% nx=2;      %number spectral elements in the model (grids)
% L=23;      %length of the medium in meters. We assumed square ->A=LxL   
% nlines=10;   %number of RCV Lines (#1 at the top, last one at the bottom)
% nrec=10;     %number of rcv for each line (#1 at the LHS)
% pml=0;  %PML boundary conditions in terms of number of spectral elements nx. If no borders write 0



dx=L/nx;
xmax=L;  zmax=L; %Model dimensions (meters)

x_pml=pml*dx; z_pml=x_pml;
x0=0+x_pml ;  xf=L-x_pml;   %first and last rcv x-positions per line (meters) 
z0=0+z_pml;  zf=L-z_pml;   %first and last rcv z-positions (meters)



%x_rcv=[x0:(xf-x0)/(nrec-1):xf];
%z_rcv=[z0:(zf-z0)/(nlines-1):zf];

%Coordinates (m)
x_rcv=[x0:(xf-x0)/(nrec-1):xf]; 
z_rcv=[z0:(zf-x0)/(nlines-1):zf];

%Matrix of the medium. 0:no station in that point   1:RCV in that point
 station=zeros(nx,nx);  
 
 %Grid number of each point
x_step=ceil(x_rcv/dx)+1;%floor(x_rcv/dx);  %Some RCV might share the same nx location. Round is not perfect. Doble check
z_step=ceil(z_rcv/dx)+1;%floor(z_rcv/dx);



%SUPER IMPORTANT- Z negative for NEW SPECFEM2D can read it the proper way!!!!
%z_rcv=-z_rcv;    %If using OLD SPECFEM@D (Anne) then, z positive (Z origin
%is in the floor of the model)

file=['STATIONS'];  %File name same as in SPECFEM2D
f=fopen(file,'w');
       
%Place the important RCV (at r,r',s,r'2) Positions by hand
formatSpec= 'S%04.0f     AA            %10.7f            %10.7f       0.0         0.0\n';
fprintf(f,formatSpec,1,8000,10000);  %RCV1 is r'
fprintf(f,formatSpec,2,6000,8000);   %RCV2 is r
fprintf(f,formatSpec,3,10000,8000);  %RCV3 is s
%fprintf(f,formatSpec,4,8000,6000);   %RCV4 is r2'(validation rcv)

rnum=3; %I placed the 1st 3 rcv "by hand" =>counter "rnum" starts at 3
  for i=1:nlines 
      for j=1:nrec
             
          
          rnum=rnum+1;   %Counting station number
                           
%Write in Specfem2D format       
formatSpec= 'S%04.0f     AA            %10.7f            %10.7f       0.0         0.0\n';
fprintf(f,formatSpec,rnum,x_rcv(j),z_rcv(i));

station(z_step(i),x_step(j))=1000;

                   end
end
            
       fclose(f); 

 
       
 %z_rcv=-z_rcv;  ONLY FOR SPECFEM V7 (New)      

%-----RCV and Source Positions for Kernel computation---------
%x starts from left to rigth.  z starts from top(surface) to bottom(matlab only)
%This is needed only if you did NOT placed the impt. rcv by hand, above 

%Position rprime in meters---
xr=8000; zr=L-10000;
%xr=4500; zr=4300.5;

 % difference (c) and location (n) in the array of my RCV position
       [c nxr] = min(abs(xr-x_rcv));
       [c nzr] = min(abs(zr-z_rcv));    
       
%Colored the area around the position to plot it      
 station(z_step(nzr-1:nzr+1),x_step(nxr-1:nxr+1))=2500;
 
 %RCV number of this position. Use it in the kernel code 
 rnum_rprime=nxr*nzr+(nrec-nxr)*(nzr-1)  %It works only if I didn't predify any rcv "by hand"
 x_step(nxr) %Position in samples
 z_step(nzr) %Position in samples
 
 
 %Position r in meters---
xr=6000; zr=L-8000;
 %xr=3000; zr=5000;

       [c nxr] = min(abs(xr-x_rcv));
       [c nzr] = min(abs(zr-z_rcv));   
 station(z_step(nzr-1:nzr+1),x_step(nxr-1:nxr+1))=2500;      
        
 %RCV number of this position. Use it in the kernel code
 rnum_r=nxr*nzr+(nrec-nxr)*(nzr-1)

 
 
       
        %Position s in meters---
xr=10000; zr=L-8000;
%xr=6000; zr=5000;

       [c nxr] = min(abs(xr-x_rcv));
       [c nzr] = min(abs(zr-z_rcv));     
 station(z_step(nzr-1:nzr+1),x_step(nxr-1:nxr+1))=2500;   
 
  %RCV number of this position. Use it in the kernel code
 rnum_s=nxr*nzr+(nrec-nxr)*(nzr-1)
 
 
 
%  
%  %----Position rprime 2 (validation RCV) in meters--
% %xr=6000-699.4755; zr=5000+1542.3;
% xr=5300.5; zr= 6542.3;
%  % difference (c) and location (n) in the array of my RCV position         
%        [c nxr] = min(abs(xr-x_rcv));
%        [c nzr] = min(abs(zr-z_rcv)); 
%        
% %Colored the area around the position to plot it      
%  station(z_step(nzr-1:nzr+1),x_step(nxr-1:nxr+1))=2900;
%  
%  
%   %----Position rprime 3 (closest to src) in meters--
% xr=5250; zr=4650.3;
% 
%  % difference (c) and location (n) in the array of my RCV position
%        [c nxr] = min(abs(xr-x_rcv));
%        [c nzr] = min(abs(zr-z_rcv));    
%        
% %Colored the area around the position to plot it      
%  station(z_step(nzr-1:nzr+1),x_step(nxr-1:nxr+1))=2500;
%  
%  
%    %----Position rprime final (far from src) in meters--
% xr=3000; zr=3601;
% 
%  % difference (c) and location (n) in the array of my RCV position
%        [c nxr] = min(abs(xr-x_rcv));
%        [c nzr] = min(abs(zr-z_rcv));    
%        
% %Colored the area around the position to plot it      
%  station(z_step(nzr-1:nzr+1),x_step(nxr-1:nxr+1))=2500;
%  
 
 

 %}


%% -------Transport & scattering mean free path computation  ---------------------------------
%{
%% Parameters - Transport & scattering mean free path computation  
nx=410;      %number spectral elements in the model (grids)
L=9000;      %length of the medium in meters. We assumed square ->A=LxL   
nlines=74;   %number of RCV Lines (#1 at the top, last one at the bottom)
nrec=1;     %number of rcv for each line (#1 at the LHS)
% We use reflecting boundaries. Thus, no PLM absorbing boundaries is needed
dx=L/nx;


x0=L/2; z0=L/2;

z_rcv=[0:L/(nlines-1):L];
x_rcv=z_rcv;



%Matrix of the medium. 0:no station in that point   1:RCV in that point
 station=zeros(nx,nx);  
x_step=ceil(x_rcv/dx)+1;%floor(x_rcv/dx);  %Some RCV might share the same nx location. Round is not perfect. Doble check
z_step=ceil(z_rcv/dx)+1;%floor(z_rcv/dx);



%SUPER IMPORTANT- Z negative so SPECFEM2D can read it the proper way!!!!
z_rcv=-z_rcv;

file=['STATIONS_cross'];  %File name same as in SPECFEM2D
f=fopen(file,'w');
       
rnum=0;

% Filling the vertical line of the cross rcv configuration
  for i=1:nlines 
     % for j=1:nrec
             
          
          rnum=rnum+1;   %Counting station number
                           
%Write in Specfem2D format       
formatSpec= 'S%04.0f     AA            %10.7f            %10.7f       0.0         0.0\n';
fprintf(f,formatSpec,rnum,x0,z_rcv(i));

%Write only coordinates
%formatSpec= '%10.7f      %10.7f\n';
%fprintf(f,formatSpec,x0,z_rcv(i));


station(z_step(i),ceil(x0/dx)+1)=1000;

                   end
%end


% Filling the Hoz line of the cross rcv configuration
 for i=1:nlines 
     % for j=1:nrec
             
          
          rnum=rnum+1;   %Counting station number
                           
%Write in Specfem2D format       
formatSpec= 'S%04.0f     AA            %10.7f            %10.7f       0.0         0.0\n';
fprintf(f,formatSpec,rnum,x_rcv(i),-z0);    %IMPORTANT -z inSPECFEM

%Write only coordinates
%formatSpec= '%10.7f      %10.7f\n';
%fprintf(f,formatSpec,x_rcv(i),-z0);


station( ceil(z0/dx)+1 ,x_step(i))=1000;  
                   end
            
       fclose(f); 

 
       
 z_rcv=-z_rcv;      
       
%}
 
 %% Visualization
       
figure(1)
 imagesc(station)
    title('RCV distribution','fontsize',16);  
    xlabel('nx','fontsize',16)
    ylabel('nz','fontsize',16)
    set(gca,'fontsize',18)
    colorbar
% [f c] = max(max(station));
 
       