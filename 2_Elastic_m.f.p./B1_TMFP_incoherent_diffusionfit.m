%Code to fit the diffusion constant D for the numerical simulations with
%SPECFEM2D 
%1 couple of coeff for each given distance 
%then statistics over the different realisations, mean, std 
clear all; close all;
clc


n=2; %dimension 2D or 3D

N_S=1;      %number source
N_R=148; %48;     %number receiver

N_C=N_S*N_R;

nb_real=40 %40; %1;   % number of realisations

load dist_cross.mat     %load distances between source-receiver
DIST=dist';

nb_dist=length(DIST);

D_test=1e6:10000:4e6; %10000:10000:1000000;      %diffusion constants to test
tau_test=1000;                   %absorption (due to the boundaries) . Now I'm not using it

D_set=zeros(nb_dist/4,nb_real);      %offsets, realisations
tau_set=zeros(nb_dist/4,nb_real); 
Delta_set=D_set; % To save missfit~errors 

I_MOY_GLOB=zeros(280,nb_dist/4,nb_real);      %temps, dist, real

dt=3e-4;
load(['wave_coh_Model1-40_Euler.mat']) %shouldn't be incoherent? \Answer: Take the envelope (coherent would be only the 1st arrival, for later times it's incoherent)
%load(['Pwave_coh_Model1-40.mat'])
%load(['Swave_coh_Model1-40.mat'])

%% AVG Data from RCV with the same offset

%----- Get Equidistant matrix -- rcv with the same offset --ADN
dist=round(dist);
%equidistant=zeros( (length(dist)-4)/4 , 4 );
equidistant=zeros( (length(dist))/4 , 4 );

for i=1:length(equidistant(:,1))  
    equidistant(i,:)=find(dist==dist(i));  %Find rcvs with the same offset
end

%plot(dist,'o')
%plot(equidistant,'o')



RESPONSE_avg_rcv=zeros( length(waves_all) , length(equidistant(:,1)) , length(waves_all(:,1,1))  );
Temp=zeros( length(waves_all) ,4 );

for real=1:40
REPONSE=squeeze(waves_all(real,:,:)); 
    
    
for i=1:length(equidistant(:,1))   
    for j=1:4
        Temp(:,j)=REPONSE(:,equidistant(i,j));  %.^2; %Matriz containing data of the equidistant rcv  Temp size is [NSTEPS, 4]
    end
    
   RESPONSE_avg_rcv(:,i, real)=mean(Temp,2);  % Avg over equidistant rcv (i.e. rcv with same offset)
 %RESPONSE_avg_rcv(:,i, real)=REPONSE(:,equidistant(i,1)); %Here, only 1 point of the 4 with same offset 

end
  

%     n=1
%     for i=1:4
%    figure (100)
%    plot( REPONSE(:,equidistant(n,i)) ) %Plot waveforms with the same offset (they look the same)
%     hold on
%     end
%     plot(RESPONSE_avg_rcv(:,n,real),'g') %Plot avg between them
    
end
    
    
%%     FIT with I_Diff

for real=1:nb_real  % realisations? (Numb Vmodels?)
    

     taille=size(RESPONSE_avg_rcv,1);  %nt
     INT=RESPONSE_avg_rcv(:,:,real).^2;  %Intensity   %RESPONSE_avg_rcv(:,:,real); 
     
    half_window=1000;           %demi fenetre moyenne glissante/half sliding window
    debut=half_window+1;    %debut du fit/beginning of the fit
    fin=taille-half_window;    %fin du fit
    pas=floor(half_window/10);   %pas du fit/not fit
    
  
    tps=[0:taille-1]'*dt; %[0:30000-1]'*dt;   %en secondes
    
    TEMPS=tps(debut:pas:fin); %time vector of sliding window
    
    I_MOY=zeros(length(TEMPS),size(INT,2));
    
    vec=debut:pas:fin;
    
    for ii=1:length(TEMPS)  %Avg slide window of Intensity
        
        I_MOY(ii,:)=mean(INT([vec(ii)-half_window:vec(ii)+half_window],:),1);     %mean intensity in the window along the time vector 
    end
    
    I_NORM=I_MOY./sum(I_MOY,1);  %I_MOY./repmat(sum(I_MOY,1),[length(TEMPS) 1]);                         %normalised intensity
    
    
    %{
    %visualization
    figure(1)  % Intensity and <I>_t
    plot(tps,INT(:,1))
    hold on
      plot(TEMPS,I_MOY(:,1),'g')
      legend('Intensity','tw Avg Intensity')
   
     % AVG Int for each offset
    for ii=1:length(I_NORM(1,:))   %size(I_NORM,2)
    
    plot(TEMPS,I_NORM(:,ii));
    title(['Avg Intensity offset ',num2str(dist(ii)),'km '])
    pause;  %Press enter to continue
    
    end
    %}
    
    
 
    
    
     DIST_GLOB=dist(1:length(equidistant));   
     DIST_GLOB=repmat(DIST_GLOB',[length(TEMPS) 1]);
     
    TEMPS_GLOB=TEMPS;
    
    DELTA=zeros( length(DIST_GLOB(1,:)),length(D_test),length(tau_test));
    
    hh=waitbar(0,'please wait...');
    
    for dd=1:length(D_test)
        for tt=1:length(tau_test)
            
            D=D_test(dd);
            tau=tau_test(tt);
            
            %I_TH=1./(4*pi*D.*TEMPS_GLOB).^(n/2).*exp(-DIST_GLOB.^2./(4*D.*TEMPS_GLOB)).*exp(-TEMPS_GLOB./tau);       %intensite theorique (equation de la diffusion), eq 5.17
            I_TH=1./(4*pi*D.*TEMPS_GLOB).^(n/2).*exp(-DIST_GLOB.^2./(4*D.*TEMPS_GLOB));      
            
            I_TH_NORM=I_TH./sum(I_TH,1);   %I_TH./repmat(sum(I_TH,1),[length(TEMPS) 1 ]);                       %intensite theorique normalisee
            
            DELTA(:,dd,tt)=sum(abs(I_TH_NORM-I_NORM),1)';                          % misfit log de valeur absolue de la difference moyenn?? sur le temps (dim1: couple; dim2: coeff D test)
                                                                                    %misfit log of absolute value of mean difference ??   on times  (dim1: couple; dim2: coeff D test)
        end
        
        waitbar(dd/length(D_test),hh);
        
    end
    
    
    close(hh);
    
    for paire=1:length(equidistant) %length(DIST)
        
        [val index]=min(squeeze(DELTA(paire,:,:)));  %minimum du misfit complet /minimum of complete misfit
        [val2 index2]=min(val);
        D_set(paire,real)=D_test(index(index2));  %Transp m.f.p.  l_star=2*D_test(index(index2))/6500
        tau_set(paire,real)=tau_test(index2);
        
        Delta_set(paire,real)=index;  %Misfit errors. ADN Oct2019
        
    end
    
    
    I_MOY_GLOB(:,:,real)=I_NORM;
    
        
end



%save resu_EL_Rx I_MOY_GLOB TEMPS D_set tau_set n DIST Delta_set

