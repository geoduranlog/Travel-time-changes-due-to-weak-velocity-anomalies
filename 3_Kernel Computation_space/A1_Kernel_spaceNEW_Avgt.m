%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             KERNEL COMPUTATION IN SPACE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear all; close all;
clc

%% ---------
%Get info from simulations   
rec(:,:,1)=load(['Several_Models/sigma5procent_NT14mil/Simulation_1/M1/S',num2str(1,'%04.0f'),'.AA.div.semc']); 
time=rec(:,1,1); %time vector (seconds)
%NSTEP=length(rec(:,1,1));  %Total time steps
dt=abs( rec(1,1,1)-rec(2,1,1) ); %dt from the simulation
 
 %---Ctes---
L=16800; %m/s
nx=865;  %410; %Nb grid points
dx=L/nx; %meters
Amax_II=1;%



%----Lame Parameters-----
vp=6500; %m/s
vs=vp/sqrt(3); %m/s
rho=3750;
mu=rho*(vs)^2;
lambda=rho*(vp^2-2*vs^2);

%lame_ratio=mu^2 / (lambda+2*mu)^2;
 Cte_p=lambda+2*mu;
  Cte_s=mu;


  %--Time Smooth parameters---
  fc=20; %Hz
  nb_period_av=8;%10;  % number of central period (1/fc) to smooth on
  half_time_av=nb_period_av*1/fc;  % duration of half the averaging window
nb_samp_half_wind=round(half_time_av/dt);

nb_sample=length(time); %I didn't do subsampling here
ind_list=nb_samp_half_wind+1:nb_sample-nb_samp_half_wind; % list of central indices for smoothing
nb_samp_smooth=length(ind_list); % number of samples of smoothed data
time_smooth=time(ind_list);


%%
%Chose Model
for M=1:4 %3 
    M
   
    
% % FOLDERS - Free Surface simulations
 folder_sim1=(['Several_Models/sigma5procent_NT14mil/Simulation_1/M',num2str(M,'%1.0f'),'/']);   %Sim1 used Exp source
 folder_sim2=(['Several_Models/sigma5procent_NT14mil/Simulation_2/M',num2str(M,'%1.0f'),'/']); % Sim2 used Exp or F45 (depends)  

 folder_save='Several_Models/sigma5procent_NT14mil/no_Lame_avgt';
 

  
%%
%====================SIMULATION I==============================
%================SOURCE AT LOCATION "S"========================  
%data recorded at r is seismogram# 2.  CTE
r=2;
i=r;
     
    r_div_simI=load([folder_sim1,'S',num2str(i,'%04.0f'),'.AA.div.semc']); 
          r_div_simI(:,2)=r_div_simI(:,2)/Amax_II; %Normalization
          
    r_curl_simI=load([folder_sim1,'S',num2str(i,'%04.0f'),'.AA.cur.semc']); 
          r_curl_simI(:,2)=r_curl_simI(:,2)/Amax_II; %Normalization
          
         
 %------P-w Intensity IP_sr-------
  IP=abs( hilbert(r_div_simI(:,2)) ).^2;
%    IP=(r_div_simI(:,2)).^2;
  

 %------S-w Intensity IP_sr-------
  IS=abs( hilbert(r_curl_simI(:,2)) ).^2;
%  IS=(r_curl_simI(:,2)).^2;
 
 %Isr=IP+IS; Isr=Isr';
 
 %Plot Raw intensities
% figure 
% plot(time,IP,'b')
% hold on
% plot(time,IS,'r')

    
%     %------Total Intensity from s to r -> u^2= ux^2+uz^2--------
   %ux=load([folder_sim1,'S',num2str(i,'%04.0f'),'.AA.BHX.semc']); 
             % ux(:,2)=ux(:,2)/Amax_II; %Normalization
              
    %uz=load([folder_sim1,'S',num2str(i,'%04.0f'),'.AA.BHZ.semc']); 
              %uz(:,2)=uz(:,2)/Amax_II; %Normalization    
         
              %u2=ux(:,2).^2 + uz(:,2).^2 ; 
              
              %---Intensities by Components----
            %  Ix=abs( hilbert(ux(:,2)) ).^2;
            %  Iz=abs( hilbert(uz(:,2)) ).^2;
    
    %Choosing only vert component
    %u2=uz;  
    
    
    %Declaration
I_SR_P_smooth=ones(nb_samp_smooth,1); I_SR_S_smooth=I_SR_P_smooth;  

    % ---time smoothing---% 
for ii=1:nb_samp_smooth
    
    ind_average=ind_list(ii)-nb_samp_half_wind:ind_list(ii)+nb_samp_half_wind;  % indices of current averaging window
    
    I_SR_P_smooth(ii,:)=mean(IP(ind_average,:),1);  % averaging
     I_SR_S_smooth(ii,:)=mean(IS(ind_average,:),1);  % averaging
     
end

%Plot Smooth Intensities
% figure 
% plot(time_smooth,I_SR_P_smooth,'b')
% hold on
% plot(time_smooth,I_SR_S_smooth,'r')
%    
  
% Declarations
nr=4903; %Total Number of rcv
Ip_save=ones(nr,length(time)); Gp_save=Ip_save;
Is_save=Ip_save;  Gs_save=Ip_save;



%% Loop over space (computing Kernel at several rprime positions)---
for i=1:nr
    
% Data recorded at r' is seismogram: arbitraty# wherever you want the kernel value

    rprime_div=load([folder_sim1,'S',num2str(i,'%04.0f'),'.AA.div.semc']);
      rprime_div(:,2)=rprime_div(:,2)/Amax_II; %Normalization
      
    rprime_curl=load([folder_sim1,'S',num2str(i,'%04.0f'),'.AA.cur.semc']); 
      rprime_curl(:,2)=rprime_curl(:,2)/Amax_II; %Normalization
    
    
%====================SIMULATION II ==============================
%================SOURCE AT LOCATION "r"========================
%data recorded at rprime is seismogram: arbitraty# wherever you want the kernel value
%i=rprime;
%Amax_II=1e10;%

    rprime_div_simII=load([folder_sim2,'S',num2str(i,'%04.0f'),'.AA.div.semc']); 
          rprime_div_simII(:,2)=rprime_div_simII(:,2)/Amax_II; %Normalization

    rprime_curl_simII=load([folder_sim2,'S',num2str(i,'%04.0f'),'.AA.cur.semc']); 
        rprime_curl_simII(:,2)=rprime_curl_simII(:,2)/Amax_II; %Normalization

    


%------P-w Intensities-------
   Ip_save(i,:)=abs( hilbert(rprime_div(:,2)) ).^2;
   Gp_save(i,:)=abs( hilbert(rprime_div_simII(:,2) )).^2;
  % Ip_save(i,:)=(rprime_div(:,2)).^2;
  % Gp_save(i,:)=(rprime_div_simII(:,2)).^2;



%------S-w Intensities-------
  Is_save(i,:)=abs( hilbert(rprime_curl(:,2)) ).^2;
  Gs_save(i,:)=abs( hilbert(rprime_curl_simII(:,2) )).^2;
%      Is_save(i,:)=(rprime_curl(:,2)).^2;
%       Gs_save(i,:)=(rprime_curl_simII(:,2)).^2;



    
end



%% =======================
Is_save=Is_save';
Ip_save=Ip_save';
Gp_save=Gp_save';
Gs_save=Gs_save';

 indicator_RawInt=(['already run Raw Intensities_M',num2str(M,'%1.0f'),])
      


%Declaration
Ip_save_smooth=ones(nb_samp_smooth,nr); Gp_save_smooth=Ip_save_smooth;
Is_save_smooth=Ip_save_smooth; Gs_save_smooth=Ip_save_smooth;
 % ---Avgt-------
 for ii=1:nb_samp_smooth
    
    ind_average=ind_list(ii)-nb_samp_half_wind:ind_list(ii)+nb_samp_half_wind;  % indices of current averaging window
    
    Ip_save_smooth(ii,:)=mean(Ip_save(ind_average,:),1);  % averaging
    Gp_save_smooth(ii,:)=mean(Gp_save(ind_average,:),1); 
  
     Is_save_smooth(ii,:)=mean(Is_save(ind_average,:),1);  % averaging
      Gs_save_smooth(ii,:)=mean(Gs_save(ind_average,:),1);
 end

%  r=2
%  figure 
% plot(time_smooth,Ip_save_smooth(:,r),'b')
% hold on
% plot(time,Ip_save(:,r),'k')
    

      
      %=====Save Raw and Smoothed Intensities=====      
 save([folder_save,'/Raw_Intensities_M',num2str(M,'%1.0f'),'.mat'],'Ip_save','Gp_save','Is_save','Gs_save','IP','IS');%,'ux','uz');  %
  save([folder_save,'/Smooth_Intensities_M',num2str(M,'%1.0f'),'.mat'],'Ip_save_smooth','Gp_save_smooth','Is_save_smooth','Gs_save_smooth','I_SR_P_smooth','I_SR_S_smooth');  %

 indicator_smoothInt=(['Saved Smooth_Intensities_M',num2str(M,'%1.0f'),])
      

  %CONVOLUTION
for i=1:nr
Nump_smooth(:,i)=conv(Cte_p.*Ip_save_smooth(:,i),Cte_p.*Gp_save_smooth(:,i)).*dt;  
Nums_smooth(:,i)=conv(Cte_s.*Is_save_smooth(:,i),Cte_s.*Gs_save_smooth(:,i)).*dt ;
end

%plot(time_smooth,Ip_save_smooth(:,1000))

% 
% %  
% r=10;
% Iptest=Ip_save_smooth(:,r:r+2);     %Iptest (nt_smooth,nr)
% Gptest=Gp_save_smooth(:,r:r+2);
% 
% for i=1:3
% Num_test1(:,i)=conv(Iptest(:,i),Gptest(:,i))  ;
% end
% 
% Num_test2=conv2(Iptest,  Gptest);     %I want conv of each column  
% Num_tes3=conv2(Gptest,1,Iptest);
  %==========================================================    
%%                  KERNEL COMPUTATION    
%==========================================================
 
%Use only 1st part of the convolution


%Kp=Nump(:,1:length(time))./(IP+IS)';
%Ks=Nums(:,1:length(time))./(IP+IS)';

%Kp=Nump_smooth(1:length(time_smooth),:)./(Cte_p.*I_SR_P_smooth+Cte_s.*I_SR_S_smooth);
%Ks=Nums_smooth(1:length(time_smooth),:)./(Cte_p.*I_SR_P_smooth+Cte_s.*I_SR_S_smooth);

Kp_hydro=Nump_smooth(1:length(time_smooth),:)./(I_SR_P_smooth); %rcv at R as Hydrophone (~acoustic =>no Lam'e parameter)
Ks_hydro=Nums_smooth(1:length(time_smooth),:)./(I_SR_P_smooth);


%
% %plot(time_smooth,Kp(:,100)-Kp_hydro(:,100))
% 
% r=3;
% figure
%   %plot(time_smooth, Nump_smooth(1:length(time_smooth),r)./(I_SR_P_smooth+0*I_SR_S_smooth))
%   %hold on 
%   plot(time_smooth, Kp(1:length(time_smooth),r),'b')
%   hold on
%     plot(time_smooth, Ks(1:length(time_smooth),r),'r')
  

 
%figure (3)
%plot(Nump(i,:),'b')

%plot(time,(IS+IP))

%========================
%%->SAVING THE KERNELS<-
%========================
folder_save='Several_Models/sigma5procent_NT14mil';
%Kernels - Expl source sim1 and Expl src sim2
  %save([folder_save,'/Kernel_M',num2str(M,'%1.0f'),'.mat'],'Kp','Ks','time_smooth');  %
  %save([folder_save,'/Kernel_hydro_M',num2str(M,'%1.0f'),'.mat'],'Kp_hydro','Ks_hydro','time_smooth');
  
 indicator_K=(['Saved Kernels_M',num2str(M,'%1.0f'),])
      


%Numerators - Expl source sim1 and Expl src sim2
%  save([folder_save,'/Nump_smooth.mat'],'Nump_smooth');  %
%  save([folder_save,'/Nums_smooth.mat'],'Nums_smooth');
% 
% 

end

toc