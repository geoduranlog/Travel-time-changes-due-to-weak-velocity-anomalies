%%% Script to determine the scattering mean free path from simulation results    (2D)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

Vp=6500;
Vs=Vp/sqrt(3)
% wave speed, determined from equipartion ratio
%v=Vp*Vs/(0.75*Vp+0.25*Vs)
v=Vp
Q=1e13;               %attenuation (large value when negligible)
nu=20;                  %central frequency

%%{
%station coordinates
folder=(['/cluster/scratch/javierd/Measure_MeanFreePaths_Elastic/RCV_CrossPattern_RX/']);

stations=importdata(['STATIONS']);
stations=stations(1,1).data;   % I take only numerical values;
stations=stations(:, [1 2]);    %Take only the (x,y) values
xr=stations(:,1);          %station coordinates x
yr=stations(:,2);          %station coordinates y 
%}
xs=8400;                 %Source position x
ys=8400;                %Source position y
nm=40;                 %number of models
nr=148;                  %number of stations
sigma=20; 

% Epicentral distances (offset of each receiver)
dist=(xr-xs).^2+(yr-ys).^2; 
dist=sqrt(dist);
x=dist;  
%save dist_cross.mat dist

load dist_cross.mat dist
x=dist;

% Some parameters
k=2*pi*nu/v; 
nstat=length(dist);

% read a file to measure tau (time step)  
ref=load(['S0001.AA.BHX.semc']);   
tau=ref(2,1)-ref(1,1);   % time step

%{
%% Method according to Frankel & Clayton, mesure maximum of envelope
% of the filtered signals and the then linear regression. 
waves_all=zeros(nm,size(ref,1),nstat);
 for Model=1:nm
    waves=zeros(size(ref,1),nstat);     %zero matrix (size ( NSTEPS, number of stations)

     for station=1:nr    %load seismogramm  
     
         %Divergence
         %s=load([folder,'M',num2str(Model),'/OUTPUT_FILES/S',num2str(station,'%04.0f'),'.AA.div.semc']);  %
         
         %Curl
         %s=load([folder,'M',num2str(Model),'/OUTPUT_FILES/S',num2str(station,'%04.0f'),'.AA.cur.semc']);  %
         
        
         %Signal=Displacement
         ux=load([folder,'M',num2str(Model),'/OUTPUT_FILES/S',num2str(station,'%04.0f'),'.AA.BHX.semc']);  %
         uz=load([folder,'M',num2str(Model),'/OUTPUT_FILES/S',num2str(station,'%04.0f'),'.AA.BHZ.semc']);  %
         s=sqrt(ux(:,2).^2+uz(:,2).^2); 
         
        
        %s=load([folder,'M',num2str(Model),'/S',num2str(35,'%04.0f'),'.AA.PRE.semp']);  
         
         
     % correction - It REMOVES geometrical and intrinsic attenuation
     s=s*sqrt(dist(station))*exp(k*dist(station)/(2*Q));    
     waves(:,station)=s;
  
     end
  waves_all(Model,:,:)=waves;
 end
 
 
%save (['wave_coh_Model1-',num2str(nm),'.mat'], 'waves_all', 'tau', 'v' ,'dist');

%}
%% Load Data
%Displacement u=abs(ux^2+uz^2)
load  (['wave_coh_Model1-',num2str(nm),'_Euler.mat'], 'waves_all', 'tau', 'v' ,'dist');  

%-------------------------------
%--OJO Use Fastest Vel so you start tw from 1st arrival---
v=Vp;
%v=Vs;
%v=Vp*Vs/(0.75*Vp+0.25*Vs);
%-------------------------------

%-Check if there are Waves with wrong values 
station=148
for k=1:nm
 find( waves_all(k,:,:)==NaN );
%a(k,:)=find(Inum(:, k, :)==NaN);
end

wave_avgM=mean(waves_all,1);

% figure 
% plot(wave_avgM(1,:,1))


figure (1)
plot(waves_all(22,:,1))
hold on
plot(waves_all(22,:,74),'--')
hold on
plot(waves_all(22,:,75),'--')
hold on
plot(waves_all(22,:,148),'--')
legend({'rcv1','rcv74','rcv75','rcv148'},'Location','northWest','FontSize',16)
%legend({'signal far offset','signal short offset'},'Location','northWest','FontSize',16)
xlim([4200,7000])
title('Waves at extreme locations same offset')

%-- Test Various Time Windows. Best are 200 or 500 samples
%TW=[20, 50, 200, 500, 1000];  
TW=500
for n=1:length(TW)
      
tw=TW(n);

maxim_all=zeros(1,nr);  maxim=maxim_all;
for Model=1:nm
    for station=1:nr
    waves=waves_all(Model,:,station);
    %waves=abs( hilbert(waves) ).^2 ;  %Use envelope if needed
    
    x=dist(station); 
    t=ceil(x/v/tau); %nt of wave 1st arrival    %ceil(x/Vp/tau)                    %beginning of first arrival (nsteps)
    
    [maxim(station),imax]=max(abs(waves(t:t+tw)));  %This if no envelope
    %[maxim(station),imax]=max(waves(t:t+tw)); %Hilbert case   %look for maximum value in a time window t:t+200, with t being the time of the P wave arrival
    end
    maxim_all=maxim_all+maxim; %Summing over Models
end
maxim=maxim_all./nm;   %The average max for each rcv


%{
figure (2)
plot( abs( hilbert(waves_all(22,:,1)) ).^2 ,'--')
hold on
plot( abs( hilbert(waves_all(22,:,74))).^2,'--')
hold on
plot(abs( hilbert(waves_all(22,:,75))).^2,'--')
hold on
plot(abs( hilbert(waves_all(22,:,148))).^2,'--')
%hold on
%plot(abs( hilbert(wave_avgM(1,:,1))).^2,'-')
%legend({'rcv1','rcv74','rcv75','rcv148','env[AvgModels_rcv1]'},'Location','northWest','FontSize',16)
legend({'rcv1','rcv74','rcv75','rcv148'},'Location','northWest','FontSize',16)
xlim([4200,7000])
title('ENV Waves at extreme locations same offset')
%}




%--Liner Regresion----
%a: slope %b: Y-axis cut
me=log(maxim);
[a,b,siga,yy,error,x]=regreslin(dist,me');                %calculate regression


% scattering mean free path obtained from regression
l(n)=-1/(2*a)  % a is slope of the regression line 
D=l*v/2 % The real Diff coef uses transp m.f. p. , not this one!
sigl(n)=abs(siga/(2*a*a))  %+-value for l
% result in Q^{-1}
%l=-2*a/k; sigl=2/k*abs(siga);  

 y=me;
% figure(3)
% plot(dist,me,'o'); hold on;
% % plot([dist(1) dist(end)],[y(1) y(end)],'b-');
% title('envelopes')


figure(4)
plot(me,'o'); hold on;
 title('Ln(max(Ampl) Vs RCVnumber)')
  xlabel('RCV Nb')



%figure
%plot(dist,'o'); hold on;

%%{
 figure (5)
  plot(dist,y,'bo',x,yy,'r-'); 
 yyp=(a+siga)*x+b;
 yym=(a-siga)*x+b;
 hold on;
 plot(x,yyp,'g-',x,yym,'k-');
 title(['Scattering m.f.p. regression'])
 xlabel('offset (m)')
 ylabel('Ln(<max|\phi(t_w)|>)')
 grid on
set(gca,'fontsize',18)
%xlim([1000,max(dist)])  %Regresion is from offset>1000m
%which grid -all

%}


end

%Mean Free time
mft=l/v