clear all
close all


load resu_EL_Rx.mat

%visu fit

position_o=1;
position_end=25;


for ii=position_o:position_end   %1:37   % ~ Loop over offsets
    
%ii=25
figure

    toto=mean(I_MOY_GLOB(:,ii,:),3); % Avg over models
    plot(TEMPS,toto);
    hold on;
    
    D=mean(D_set(ii,:));   % Avg over models  %%Transp m.f.p.  l_star=2*D/6500
    tau=mean(tau_set(ii,:));
    
    
    I_TH=1./(4*pi*D.*TEMPS).^(n/2).*exp(-DIST(ii).^2./(4*D.*TEMPS)); %.*exp(-TEMPS./tau);       %intensite theorique (equation de la diffusion)
    
    I_TH_NORM=I_TH./sum(I_TH);                       %intensite theorique normalisee
    
    plot(TEMPS,I_TH_NORM,'r');
 title(['Intensities at offset ' num2str(round(DIST(ii))),'km'])
legend({'Avg Intensity','Itheo'},'Location','northWest','FontSize',16) 
grid on 
xlabel('time(s)')
set(gca,'fontsize',18)
    
    pause;
    %hold off
    
end

    

x=DIST(ii) % Offset of that position
%x_max=(DIST(position_o))    %Max offset used (m)
%x_min=(DIST(position_end))    %Min offset used (m)



%--Use Diff Coeff only from rcv with a 'proper' offset --
%  Too short offset, there is no diffusion
% Too far, there migth be boundaries effects
D_good=D_set(1:25,:); % Choose 'proper offsets'
D_goodModels=mean(D_good,2);
D_good=mean(mean(D_good,2)); %%Avg over Models and offsets  %2*D_good/6500  

%TEMP
%DIST_good=DIST(position_o:position_end);

 
Delta_good=mean( mean(Delta_set(1:25,:),2) ); %Error of Diff Coef  ADN Oct2019

Vp=6500;
Vs=Vp/sqrt(3)
% wave speed, determined from equipartion ratio
c=Vp*Vs/(0.75*Vp+0.25*Vs)  %Using E velocity

%%======Transp m.f.p.====== 
l_star=2*D_good/c  %l_star=2*D_goodModels/c     save Diff_withDist D_goodModels DIST_good

l_error=(2/c)*Delta_good   

l_star*l_error  %Error