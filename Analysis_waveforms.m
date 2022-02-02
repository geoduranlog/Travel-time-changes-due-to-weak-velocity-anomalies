%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             EXTRACT INFORMATION FROM SPECFEM2D OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%-----Read geophone locations (Stations)-----
%X location is in the first Col, i.e, stations(:,1)
stations=importdata('OUTPUT_FILES/STATIONS');
%To acces the variables inside the struct: structName(indices).fieldName
stations=stations(1,1).data;   % I take only numerical values;
stations=stations(:, [1 2]);    %Take only the (x,y) values

%----NORMALIZATION Parameter!!-----
%Normalization to simulate an elastic unit force - If you don't need it
%just set it as 1
Amax_II=1; %1.7649e-13; %1.7419e-13;     %% i.e. Amax of SimulationII rcv1 ->  m(1)


%-----Read all seismograms----------
%%---Vertical Component--
for i=1:length(stations(:,1))  %i=1:Nb geophones
          
              %num2str(i,'%04.0f') gives me the apropiate values when reading any seismogram                         
      rec(:,:,i)=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.BXZ.semc']); 
      
end

%---Reading curl(displacement)--- 
for i=1:length(stations(:,1))  %i=1:Nb geophones

rec_curl(:,:,i)=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.cur.semc']); 
      
end

%---Reading Divergence of displacement--- 
for i=1:length(stations(:,1))  %i=1:Nb geophones
    
            rec_div(:,:,i)=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.div.semc']); 
      
end


%%--Horizontal Component---
% for i=1:length(stations(:,1))
%        rec(:,:,i)=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.BXX.semc']); 
%  end



xstations=stations(:,1);   %Take only x-position of the stations
for j=1:length(xstations)    %There are 10 stations
for i=1:length(rec(:,1,1))  %Wavelet, There are 1500 values of time and amplitude.
x(i,j)=xstations(j)';  %I save the position of the station, for each station.
end
end



%---Save the Max_ampl of each recording signal and its value of time in s.
%Check which component you took few lines above, Hoz or Vertical.
for i=1:length(stations(:,1))  %i=1:Nb geophones
  [m(i), n(i)] = max( abs(rec(:,2,i)) );   % [maxAmp_value  n_position_in_t_vector]
  n(i)=rec(n(i),1,i);               %  n_position_in_t_vector -> I save the t_value
end


% 
% %---Save the first non zero Amp value of each recording signal and also its value of time in s.
% for i=1:length(stations(:,1))  %i=1:Nb geophones
%     
%     %It doesn't work if the geophone doesn't record at least one a non zero
%     %Amp. value ->Fix it later
%     
%     k(i,1)=find(rec(:,2,i)~=abs(0),1); %Find position of 1st nonzero Amp. 
%     k(i,2)=rec(k(i),2,i);        %Save  its value of Ampl       
%     k(i,1)=rec(k(i),1,i);        %Save its value of time  => Finally  % k (nposition_first_nonzero_Amp, value_of_its_time)  
% end




%------Histogram of Min Nb. of points per wavelenght---------
%Conditions to be satisfied:
%For ACUSTIC Media --> Nb. of points>= 5.5
%For ELASTIC Media --> Nb. of points>= 4.5

figure(1)
hist=load('OUTPUT_FILES/points_per_wavelength_histogram_S_in_solid.txt');
bar(hist(:,1),hist(:,2))
grid on
title('Points per wavelenght','fontsize',14);
xlabel('Range of min number of points per S wavelength in solid')
ylabel('Percentage of elements (%)');

% 
% figure(2)
% hist=load('OUTPUT_FILES/points_per_wavelength_histogram_P_in_fluid.txt');
% bar(hist(:,1),hist(:,2))
% grid on
% title('Points per wavelenght','fontsize',14);
% xlabel('Range of min number of points per P wavelength')
% ylabel('Percentage of elements (%)');

%% -----------FIGURES-----------------------
%%--Plot all seismograms
% figure(3)
% for i=1:length(stations(:,1))
% plot3(x(:,i),rec(:,1,i),rec(:,2,i))
% hold on
% end
% hold on
% %plot3(xstations,n,m,'k-')  %Plot Recorded Signals Max Amp
% hold on
% %plot3(stations,k(:,1),k(:,2),'k-')  %Plot Recorded Signals 1st non zero
% grid on
% title('Recorded Signals','fontsize',14);
% xlabel('X (m)')
% ylabel('Time (s)');
% zlabel('Amplitude (m)');

% %--Plot all "Curl" seismograms
% figure(23)
% for i=1:length(stations(:,1))
% plot3(x(:,i),rec_curl(:,1,i),rec_curl(:,2,i))
% hold on
% end
% hold on
% %plot3(xstations,n,m,'k-')  %Plot Recorded Signals Max Amp
% hold on
% %plot3(stations,k(:,1),k(:,2),'k-')  %Plot Recorded Signals 1st non zero
% grid on
% title('Recorded Curl(u) Signals','fontsize',14);
% xlabel('X (m)')
% ylabel('Time (s)');
% zlabel('Amplitude (m)');
% 
% %--Plot all "Div" seismograms
% figure(33)
% for i=1:length(stations(:,1))
% plot3(x(:,i),rec_div(:,1,i),rec_div(:,2,i))
% hold on
% end
% hold on
% %plot3(xstations,n,m,'k-')  %Plot Recorded Signals Max Amp
% hold on
% %plot3(stations,k(:,1),k(:,2),'k-')  %Plot Recorded Signals 1st non zero
% grid on
% title('Recorded Div.(u) Signals','fontsize',14);
% xlabel('X (m)')
% ylabel('Time (s)');
% zlabel('Amplitude (m)');


%---------Analysis 1st Receiver--------
first_rcv=1;
i=first_rcv;
    %Check  specific data in Results_SpecFem2D
    
datax=load(['OUTPUT_FILES/AA.S000',num2str(i),'.BXX.semc']);

dataz=load(['OUTPUT_FILES/AA.S000',num2str(i),'.BXZ.semc']);
%----NORMALIZATION-----
%Normalization to simulate an elastic unit force - Green's function
%calculation in elastic media
datax(:,2)=datax(:,2)/Amax_II;   dataz(:,2)=dataz(:,2)/Amax_II;
%max(max(abs(dataz(:,2))))


%----Plot Hoz and Vert Component- First Receiver
figure(15) 
subplot(2,1,1) 
plot(datax(:,1), datax(:,2))
%xlim([-2 10]);
title(['Signal Receiver ' num2str(first_rcv), ' Horizontal Comp.'],'fontsize',14); grid on
%title('Signal First Receiver Horizontal Comp.','fontsize',14);grid on
subplot(2,1,2) 
plot(dataz(:,1), dataz(:,2))
title(['Signal Receiver ' num2str(first_rcv), ' Vertical Comp.'],'fontsize',14);
%title('Signal First Receiver Vertical Comp.','fontsize',14);
grid on
%xlim([-2 10]);




%---------Analysis Last Receiver--------
 %i=length(n);
 last_rcv=2;
  i=last_rcv;
          datax=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.BXX.semc']);
          dataz=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.BXZ.semc']);        
%----NORMALIZATION-----
%Normalization to simulate an elastic unit force - Green's function
%calculation in elastic media 
datax(:,2)=datax(:,2)/Amax_II;   dataz(:,2)=dataz(:,2)/Amax_II;
          
%----Plot Hoz and Vert Component- Last Receiver
figure(16) 
subplot(2,1,1) 
plot(datax(:,1), datax(:,2)) 
%title('Signal Last Receiver Horizontal  Comp.','fontsize',14); grid on
title(['Signal Receiver ' num2str(last_rcv), ' Horizontal Comp.'],'fontsize',14); grid on
subplot(2,1,2) 
plot(dataz(:,1), dataz(:,2))
%title('Signal Last Receiver Vertical  Comp.','fontsize',14);
title(['Signal Receiver ' num2str(last_rcv), ' Vert. Comp.'],'fontsize',14);
grid on
 
m(i)



%---
 m(i)  %-> Max AMplitude value for station 'i' (Check if it comes from
% Vertical or Hz component!! -> first part of this code)
% 
% mean(dataz(:,2))
% 
% proc=load(['proc000000_rho_vp_vs.dat']);
% Vp=proc(:,4);
% %find(Vp~=a,1)-> Find the position where there is the first value which is different to 'a' 
% find(Vp~=6382.5,1)  %6273.7
% 
% Vs=proc(:,5);
% find(Vs~=4021.5,1)

%%
%===========Analysis Curl(u)==============

%-----First Receiver------------
i=first_rcv;
    %Check  specific data in Results_SpecFem2D
    data_curl=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.cur.semc']);
        data_curl(:,2)=data_curl(:,2)/Amax_II;  %Normalization

%----Plot Curl - First Receiver
figure(25) 
%subplot(2,1,1) 
plot(data_curl(:,1), data_curl(:,2))
title(['Signal Receiver ' num2str(first_rcv), ' Curl(u)'],'fontsize',14);
%title('Signal First Receiver Curl(u)','fontsize',14);grid on
grid on

clear data_curl

%-----Last Receiver------------
i=last_rcv;
    %Check  specific data in Results_SpecFem2D   
 data_curl=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.cur.semc']);
  data_curl(:,2)=data_curl(:,2)/Amax_II;  %Normalization

%----Plot Curl - Last Receiver
figure(26) 
%subplot(2,1,1) 
plot(data_curl(:,1), data_curl(:,2))
title(['Signal Receiver ' num2str(last_rcv), ' Curl(u)'],'fontsize',14);
%title('Signal Last Receiver Curl(u)','fontsize',14);grid on
grid on

%%

%===========Analysis Div.(u)==============

%-----First Receiver------------
i=first_rcv;
    %Check  specific data in Results_SpecFem2D
         data_div=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.div.semc']);
          data_div(:,2)=data_div(:,2)/Amax_II;  %Normalization
 
          data_div_rcv1=data_div;
          
%----Plot Div - First Receiver
figure(35) 
%subplot(2,1,1) 
plot(data_div(:,1), data_div(:,2))
title(['Signal Receiver ' num2str(first_rcv), ' Div.(u)'],'fontsize',14);
%title('Signal First Receiver Div.(u)','fontsize',14);grid on
grid on

%clear data_div
%-----Last Receiver------------
i=last_rcv;
    %Check  specific data in Results_SpecFem2D
            data_div=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.div.semc']);
             data_div(:,2)=data_div(:,2)/Amax_II;  %Normalization

%----Plot Div - Last Receiver
figure(36) 
%subplot(2,1,1) 
plot(data_div(:,1), data_div(:,2))
title(['Signal Receiver ' num2str(i), ' Div.(u)'],'fontsize',14);
%title('Signal Last Receiver Div.(u)','fontsize',14);grid on
grid on




%------Total Intensity from s to r -> u^2= ux^2+uz^2--------
    ux=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.BXX.semc']); 
              ux(:,2)=ux(:,2)/Amax_II; %Normalization
              
    uz=load(['OUTPUT_FILES/AA.S',num2str(i,'%04.0f'),'.BXZ.semc']); 
              uz(:,2)=uz(:,2)/Amax_II; %Normalization    
              
        u2=ux(:,2).^2 + uz(:,2).^2 ; 
        
        u=sqrt(u2);
   
         method='linear';
[upperenv lowerenv] = envelope(u2, method);    
I_stor=upperenv; 
        
        
        ratio=data_div_rcv1(:,2)./u;
        
        ratio2=(data_div_rcv1(:,2).^2)./(u.^2);
        
 
        
        nto=1; %40224;
figure(37);
%plot(data_div(nto:end,1), data_div_rcv1(nto:end,2).^2)
%hold on
%plot(data_div(nto:end,1), u(nto:end))
%hold on
plot(data_div(nto:end,1), ratio(nto:end))
hold on
plot(data_div(nto:end,1), ratio2(nto:end))
title(['Large travel time behaviour'],'fontsize',14);
    %legend({'\nabla.u ->Intensity s to rprime ','u -> I total'},'FontSize',14)
    legend({'\nabla.u / u', 'ratio Isr"/ Itotal'},'FontSize',14) 
    xlabel('Travel Time (s)')
   grid on

   
   figure(38);
%plot(u(nto:end,1),dataz(:,2))
plot(dataz(:,1),u)
title(['Displacement vs time'],'fontsize',14);
    %legend({'\nabla.u ->Intensity s to rprime ','u -> I total'},'FontSize',14)
    legend({'u '},'FontSize',14) 
    xlabel('Travel Time (s)')
    ylabel('Displacement (units)')
   grid on
   
     figure(39);
plot(dataz(:,1),I_stor)
title(['Total Intensity'],'fontsize',14);
    %legend({'\nabla.u ->Intensity s to rprime ','u -> I total'},'FontSize',14)
    legend({'I_{T} '},'FontSize',14) 
    xlabel('Travel Time (s)')
    ylabel('Total Intensity')
   grid on