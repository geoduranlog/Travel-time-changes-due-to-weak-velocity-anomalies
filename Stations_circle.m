%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      CREATE STATIONS (FILE) CIRCULAR PATTERN FOR SPECFEM2D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Useful when you want to plant many receivers along a circle around the
%perturbation point

close all
clear all

%Parameters
R=1000; %2000;  %Radius (m) of circunference 
vp=6500;  vs=vp/sqrt(3);
lambda=vp/20; % P-wavelength . Source fo=20Hz
lambda_s=vs/20; % S-wavelength . Source fo=20Hz

%Anomaly coordinates (m)
%xo=8400.24; yo=6300.18;     %   6533.52    6300.18    3750.00 
xo=2000.00; yo=2000.00;   

%%
%a=4*pi*R/lambda;  % Factor by which I divide the circunference when I want to obtain arc=R*dTheta~lambda/2

%dtheta=2*pi/a;
dtheta=lambda_s/(3*R);   %Place min 2 RCV per min. lambda (i.e. 1rcv every lambda/2)
theta=0:dtheta:2*pi; %2*pi-dtheta;  %Angle vector in Radians, Angle 0 and direcction follows usual Cartesian coordinates


%RCV Coordinates
x=xo+R*cos(theta);
y=yo+R*sin(theta);

%% -----RCVs Fill the circle ----
%%{

%Circle - coordinates
%file=['STATIONS_circle'];  %File name same as in SPECFEM2D
%file=['STATIONS_circle_lambda_S']; 
f=fopen(file,'w');
       
rnum=0;

  for i=1:length(x) %nlines 
      %for j=1:nrec
             
          
          rnum=rnum+1;   %Counting station number
                           
%Write in Specfem2D format       
formatSpec= 'S%04.0f     AA            %10.7f            %10.7f       0.0         0.0\n';
fprintf(f,formatSpec,rnum,x(i),y(i));

%station(z_step(i),x_step(j))=1000;

                   %end
end
            



%% 
%{
% Use as a reference a previous Vel model
%it contains already the correct dimensions of the medium
M=1;
cd (['/Users/Alejandro/KernelsComparison_paper/Normalization_factor/OldSpecfem/Hetero_Elastic/Models/M',num2str(M,'%01.0f')]);


%Read ASCII-delimited file of numeric data into matrix
Vmodel = dlmread('VmodelEl_sigma20_Vp_pert6955_r5.xyz');
cd ../../

%Line of that point in the text file
line=find(   round(Vmodel(:,1)==xo) &  round(Vmodel(:,2)==yo)  );


d=round(lambda/2);

    %Possible x-coord of the circle 
for n=0:100
    if xo-R+n*d<=xo+R
x(n+1)=xo-R+n*d;

    else
        break
    end
end

y1= sqrt(R^2 - (x-xo).^2) +yo; %Positive semicircle
y2= -sqrt(R^2 - (x-xo).^2) +yo; %Negative semicircle



%% -----RCVs Fill the circle ----
%%{

%Positive semicircle - coordinates
file=['STATIONS_circle'];  %File name same as in SPECFEM2D
f=fopen(file,'w');
       
rnum=0;

  for i=1:length(x) %nlines 
      %for j=1:nrec
             
          
          rnum=rnum+1;   %Counting station number
                           
%Write in Specfem2D format       
formatSpec= 'S%04.0f     AA            %10.7f            %10.7f       0.0         0.0\n';
fprintf(f,formatSpec,rnum,x(i),y1(i));

%station(z_step(i),x_step(j))=1000;

                   %end
end
            
    
       
       %Negative semicircle - coordinates

  for i=1:length(x) %nlines 
      %for j=1:nrec
             
          
          rnum=rnum+1;   %Counting station number
                           
%Write in Specfem2D format       
formatSpec= 'S%04.0f     AA            %10.7f            %10.7f       0.0         0.0\n';
fprintf(f,formatSpec,rnum,x(i),y2(i));

%station(z_step(i),x_step(j))=1000;

                   %end
end
            
       fclose(f); 
       
       
%} 
       

