
clear all
close all


%Get info from simulations 
folder_sim1='Several_Models/sigma5procent_NT14mil/Simulation_1/M1/'; % _Exp  _Dirac
 rec(:,:,1)=load([folder_sim1,'S',num2str(1,'%04.0f'),'.AA.div.semc']);    
time=rec(:,1,1); %time vector (seconds)
%NSTEP=length(rec(:,1,1));  %Total time steps
 dt=abs( rec(1,1,1)-rec(2,1,1) ); %dt from the simulation

 nr=4903; %total rcv
 nt=14000; %6500;  %total time steps (smooth_time)
 nm=1 %4%7;  %total models
 
 
 %----Lame Parameters-----
vp=6500; %m/s
vs=vp/sqrt(3); %m/s
rho=3750;
mu=rho*(vs)^2;
lambda=rho*(vp^2-2*vs^2);

%lame_ratio=mu^2 / (lambda+2*mu)^2;
 Cte_p=lambda+2*mu;
  Cte_s=mu;
  
  %Load only time_smooth
  %load(['Several_Models/sigma5procent_NT14mil/Kernel_M1_old']); clear Kp Ks
   %--Time Smooth parameters---
  fc=20; %Hz
  nb_period_av=8;%10;  % number of central period (1/fc) to smooth on
  half_time_av=nb_period_av*1/fc;  % duration of half the averaging window
nb_samp_half_wind=round(half_time_av/dt);

nb_sample=length(time); %I didn't do subsampling here
ind_list=nb_samp_half_wind+1:nb_sample-nb_samp_half_wind; % list of central indices for smoothing
nb_samp_smooth=length(ind_list); % number of samples of smoothed data
time_smooth=time(ind_list);
 %--------

 %%Smooth Intensities
%IP=ones(nt,nr,nm); IS=IP; GP=IP;  GS=IP;  Isr_P=ones(nt,nm); Isr_S=Isr_P;

 %Raw Intensities
Ip=ones(nt,nr,nm); Is=Ip; GP=Ip;  GS=Ip;  Isr_P=ones(nt,nm); Isr_S=Isr_P; 

%--Load Kernels in space (no norm)- All models--
for M=1:nm
    
    %--Case: Smooth Intensities ---
%load(['Several_Models/no_Lame_avgt/Smooth_Intensities_M',num2str(M,'%1.0f')]); 

% IP(:,:,M)=Ip_save_smooth;  
% GP(:,:,M)=Gp_save_smooth;
% IS(:,:,M)=Is_save_smooth;    
% GS(:,:,M)=Gs_save_smooth;  
% 
% Isr_P(:,M)=I_SR_P_smooth;
% Isr_S(:,M)=I_SR_S_smooth;


%--Case: Raw Intensities ---
load(['Several_Models/sigma5procent_NT14mil/no_Lame_avgt/Raw_Intensities_M',num2str(M,'%1.0f')]); 
 
Ip(:,:,M)=Ip_save;  
GP(:,:,M)=Gp_save;
Is(:,:,M)=Is_save;    
GS(:,:,M)=Gs_save;  

Isr_P(:,M)=IP;     %The notation here is different form the Smoothed case
Isr_S(:,M)=IS;

%clear Ip_save_smooth Gp_save_smooth ...

end
clear IP IS  
IP=Ip; IS=Is;  %For Raw intensities only


% %Skip some models
% mo=3; mf=6; %models which remain
% IPskip= IP(:,:,mo:mf);
% GPskip= GP(:,:,mo:mf);
% ISskip= IS(:,:,mo:mf);
% GSskip= GS(:,:,mo:mf);
% 
% Isr_Pskip=Isr_P(:,mo:mf);
% Isr_Sskip=Isr_S(:,mo:mf);


%Avg over models
IPmean=mean(IP,3);
GPmean=mean(GP,3);
ISmean=mean(IS,3);
GSmean=mean(GS,3);

Isr_Pmean=mean(Isr_P,2);
Isr_Smean=mean(Isr_S,2);

% IPmean2=mean(IPskip(:,:,3:4),3);
% GPmean2=mean(GPskip(:,:,3:4),3);
% ISmean2=mean(ISskip(:,:,3:4),3);
% GSmean2=mean(GSskip(:,:,3:4),3);
% 
% Isr_Pmean2=mean(Isr_Pskip(:,3:4),2);
% Isr_Smean2=mean(Isr_Sskip(:,3:4),2);



%==Save Imean==
% save(['Several_Models/sigma5procent_NT14mil/no_Lame_avgt/Mean_Intensities.mat'],'IPmean','GPmean','ISmean','GSmean','Isr_Pmean','Isr_Smean');  %
%load (['Several_Models/sigma5procent_NT14mil/no_Lame_avgt/Mean_Intensities.mat'],'IPmean','GPmean','ISmean','GSmean','Isr_Pmean','Isr_Smean');  %


  %CONVOLUTION
for i=1:nr
Nump(:,i)=conv(Cte_p.*IPmean(:,i),Cte_p.*GPmean(:,i)).*dt;  
Nums(:,i)=conv(Cte_s.*ISmean(:,i),Cte_s.*GSmean(:,i)).*dt ;
end



time_smooth=time;

a=0;  %include(a=1) or remove(a=0) "Is" in the denominator of the kernel
%Kp=Nump(1:length(time_smooth),:)./(Cte_p.*Isr_Pmean+a*Cte_s.*Isr_Smean);
%Ks=Nums(1:length(time_smooth),:)./(Cte_p.*Isr_Pmean+a*Cte_s.*Isr_Smean);

 %Convention->you recorded with an Hydrophone in R (~acoustic => no Lame parameter in the denominator)
Kp=Nump(1:length(time_smooth),:)./(Isr_Pmean); 
Ks=Nums(1:length(time_smooth),:)./(Isr_Pmean);



 
% plot(Nums(:,1))
 
figure
plot(time_smooth,Nump(1:length(time_smooth),2)./(Cte_p.*Isr_Pmean+0*Cte_s.*Isr_Smean) )

%% ---   Image  ---
x0=0; z0=x0;
nrec=70;  nlines=70;
nx=865;
L=16800;
dx=L/nx;
xf=L; zf=xf; %m

x_rcv=[x0:(xf-x0)/(nrec-1):xf];
z_rcv=[z0:(zf-x0)/(nlines-1):zf];


%Matrix of the medium. 
 kernelp=zeros(nlines,nrec);   kernels=kernelp;
%x_step=ceil(x_rcv/dx)+1;%floor(x_rcv/dx);  %Some RCV might share the same nx location. Round is not perfect. Doble check
%z_step=ceil(z_rcv/dx)+1;%floor(z_rcv/dx);

%--Adapt--
time=time_smooth;



%============================
%Choose travel time<time(end)
%============================
travelt=1.75 %0.75 %seconds

travelt=round( (travelt-time(1))/dt +1); %samples

if travelt>length(time)
    travelt=length(time);
end

if travelt<1
    travelt=1;
end
 
 %%
 
 
 

 
count=4;   %Start from 4 because first 3rcv will be placed "by hand"
for j=1:70
    for i=1:70 

     
         
   if count<=nrec*nlines+3   %Total number of stations
        
        
    kernelp(j,i)=Kp(travelt,count) ;
    kernels(j,i)=Ks(travelt,count); 
    
  %kernelp(j+(70-j),i)=Kp(travelt,count) ;
  %kernels(j+(70-j),i)=Ks(travelt,count); 
   
  % kernelp(70-(70-j),i)=Kp(travelt,count) ;
   %kernels(70-(70-j),i)=Ks(travelt,count); 
    
   
   count=count+1;
   end
   
   end
    
    
end

     %Free surface is in the top of the image 
     kernelp=flipud(kernelp) ;   % kernelp=flip(kernelp) ;
     kernels=flipud(kernels);   % kernels=flip(kernels);
 
     
       
       %% NORMALIZATION FACTOR
% It is the injected Energy of the impulse
%load /Users/alejandro/Kernel_ElasticMedia_Project_Specfem/Matlab_Codes_Elastic/Kernel_in_space/norm_p_sim2_srcE 
% load /Users/alejandro/Kernel_ElasticMedia_Project_Specfem/Matlab_Codes_Elastic/Kernel_in_space/norm_s_sim2_srcE
%norm_green_p=max(norm_green_p);
%norm_green_s=max(norm_green_s);
%norm_p=max(norm_p);
%norm_s=max(norm_s);
%NORMFACTOR= norm_green_s+norm_green_p;
NORMFACTOR=1e20;  %2.9921e+18; %1; %1e20=max(I(t) Explsrc sim2 at position r)= Afactor^2;

kernelp=kernelp./NORMFACTOR;

kernels=kernels./NORMFACTOR;
       
       

    %Interpolation
k_dense=interp2(kernels,6);
%k_dense=kernels;
  
%%
%===Figures========

figure(1)
% imagesc(x_rcv, z_rcv ,k_dense)
 surf(z_rcv/1000, x_rcv/1000 , kernels)
    title('Kernel KS in space','fontsize',16);  
    xlabel('x (km)','fontsize',16) %nx
    ylabel('z (km)','fontsize',16)
    set(gca,'fontsize',18)
    colorbar
     colormap(autumn(5)) %autumn(5))

  

%%


% folder_sim1='Several_Models/sigma5procent_NT14mil/Simulation_2/M1/'; 
% 
% u_R=load([folder_sim1,'S',num2str(2,'%04.0f'),'.AA.cur.semc']);
% u_S=load([folder_sim1,'S',num2str(3,'%04.0f'),'.AA.cur.semc']);
% 
% IS_R=abs( hilbert(u_R(:,2)) ).^2;
% IS_S=abs( hilbert(u_S(:,2)) ).^2;
% 
% 
% 
% %figure(3)
% %plot(time,up_R(:,2),'r')
% %hold on
% %plot(time,up_S(:,2),'k')
% %ylim([-1e-6,1e-6])
% %ylim([-15e-5,7e-5]) %For Ks  
% 
% figure(4)
% plot(time,IS_R,'r')
% hold on
% plot(time,IS_S,'k')
% 
% 
% max(up_S(:,2))