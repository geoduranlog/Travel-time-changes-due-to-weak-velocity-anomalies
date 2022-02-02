% random 2D velocity model after Frankel & Clayton (1986)
clear all; 
close all;
tic
%%%% Model description %%%%
% size of the medium   % if you change the size, change also nspec
xmax=33600; %16800 ; %in m
ymax=33600; %16800 ; %in m  

% Properties of the medium 
correlation='von_karman'; %'gauss' or 'von_karman'
ax=300; ay=300; % correlation length
sigma=0.2;%      %0.050 -> velocity perturbations = +- 5%
rho=3750;  %2700;  %%kg/m3

% Mesh size
% Elle doit ??tre un peu plus grande que le maillage de specfem2d
dx=20; dy=20;  % dx=dz=19.5 m


% GLL points
nspec=432; %216;  %54;
xgll0=[0 13.43 38.89 64.35 77.78]; %pour Mnx=54, max=4200 !!
% ?? lire dans le fichier gll.txt 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other parameters
% dans l'espace distance
xbuff=2*ax; % we enlarge the space to avoid Gibbs -- ugly ?
x=[-xbuff:dx:xmax+xbuff];
nx=length(x); 

ybuff=2*ay; 
y=[-ybuff:dy:ymax+ybuff];
ny=length(y);

% dans l'espace nombre d'onde  
dkx=1/((nx-1)*dx).*2*pi;
kxmax=1/dx.*2*pi;
kx=[-kxmax/2:dkx:kxmax/2];

dky=1/((ny-1)*dy).*2*pi;
kymax=1/dy.*2*pi;
ky=[-kymax/2:dky:kymax/2];



%------------------------------------------------------------------------
%----------------------LOOP FOR SEVERAL MODELS----------------------------
%---How many Vel Models----
        nmodels=60%20% 35% 40    
        
        %vs_min=zeros(1,nmodels); 
         count=0; % Count if there are vs<0 

        for M=54:60%2:nmodels %20 %30 %20%1:10 %35:nmodels
       count=count+1;
            
%%Go to the model directory 
cd (['Models/M',num2str(M,'%01.0f')]);

%%{
%%%% Construction of the random model in space 
z=random('norm',0,sigma,ny,nx); % loi normale autour de zero

% TF 2D, shift et spectre de puissance
Z=fftshift(fft2(z));  
Pzz=abs(Z.*conj(Z));


% Pour un milieu avec ax neq ay, kr^2 = kx^2 + ky^2

switch correlation
    case {'gauss'}
     % Construction du filtre : Gauss 
     for i=1:nx 
        for j=1:ny
            Pfil(j,i)=ax*exp(-kx(i)*kx(i)*ax*ax/4)*...
                ay*exp(-ky(j)*ky(j)*ay*ay/4)/2;
        end
     end
      
    case {'von_karman'}
    % Construction du filtre : Von Karman
    for i=1:nx
        for j=1:ny
            Pfil(j,i)=ax*ay/(1+kx(i)*kx(i)*ax*ax+ky(j)*ky(j)*ay*ay);
        end
    end
end


%%% Application du filtre
Pzz_fil=Pzz.*Pfil;

%%% Retour en x,y  (ne pas oublier le ifftshift, on suppose phase=0)
z_fil=real(ifft2(sqrt(ifftshift(Pzz_fil))));
%}

%%% Projection du mod??le sur les GLL points
xgll=xgll0; ngll=length(xgll);
for ispec=2:nspec
    xgll=[xgll(1:ngll) xgll0(2:5)+xgll(ngll)];
    ngll=length(xgll);
end
ygll=xgll;
[Xgll,Ygll] = meshgrid(xgll,ygll);
[X,Y] = meshgrid(x,y);


% Coment Filter construction
%%{  
zgll = interp2(X,Y,z_fil,Xgll,Ygll,'linear'); 
% linear donne plus de ressemblance (visuellement) entre 
% z_fil et z_gll


%%% Manipulation du signal  -- LAID ?
moy_target=0;
et_target=sigma;

% Chgt de la moyenne
moy=mean(zgll(:));    imagesc(zgll); colorbar
z_final=zgll+moy_target-moy;    imagesc(z_final); colorbar

% Chgt de l'??cart type  %%Coment this part ONLY if you want to read the
% non-pert model and place a localize anomaly as it is showed few lines below
et=std(z_final(:));
z_final=z_final.*et_target/et; %%% z = dv/v0


% OJO All elemets should be >-1, so I have positive Velocity when doing vref*(1+z_final). (EL media). ADN January 2020
idx=find(z_final<-1);
z_final(idx)=(-0.9-0.9).*rand(size(z_final(idx)))+0.9;  % Replace by random num within [-0.9 0.9] range



 %}  %End Coment Filter construction
%========================
%========================
%->SAVING THE MODEL<-
 %save(['VmodelLarge_M',num2str(M,'%01.0f'),'.mat'],'z_final');

%========================
%========================
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%% ---- LOAD NON PERTURBED MODEL ------
%clear z_final
%load(['VmodelLarge_M',num2str(M,'%01.0f'),'.mat']);   

%% LAYER PERTURBATION
% write velocities on the model/ possibility to include perturbed layers
%%{
         vref=6500;   %background velocity
         v1=vref+100000; %vref*(1 + 7/100) ;     %perturbed velocity
         z_final1=vref*(1+z_final);    %background
         z_final2=v1*(1+z_final);      % perturbed layer
         
         % if z_final1<0 || z_final2<0
         %if NON empty
           %  [a,b]=findz_final1;   z_final1(a)=z_final1(a+1) end ...
  
         %---Depth of perturbed layer---
         %Specfem depth location (Convention Specfem is z starts at the base of
         %the medium)
         %de_Spec=3953.35; %6300.18  %meters
         %depth=ymax-de_Spec;  %Change to Matlab convention (z starts at the top)
       depth=11000;
         
         thickness=0;  %thickness of perturbed layer (m)  Write 0 if no layer anomaly
         nygll1=min(find(ygll>=depth))-2;   % Test it ygll(nygll1) ~ depth
        y_bottom=ygll(nygll1)+thickness;  %Bottom in Matlab convention
        [minValue,nygll2] = min(abs(y_bottom-ygll)) % Test it ygll(nygll2)=depth+thickness
         
                 %nygll1=round(depth/ymax*ngll);  %217=4200  868=16800
                 %nygll2=round((depth+thickness)/ymax*ngll);
         
                 
         
          
%}

%{
                   %============Small Perturbation=========
  %--Position in x    
        lx=min(find(xgll>=869));  nxo=lx;      %xgll(lx) %r5_x=869.01
        %nxo=round(length(z_final1)/2); 
 
 %--Position in z (far field i.e. z>=4*trans. mean free path=4*940
  %near field ~trans.m.f.p./2)
        lz=min(find(ygll>=6300));  nzo=lz;  %ygll(lz)   r5_z=6300.18 
        %nzo=round(length(z_final1)/2); 
 
   c=0;  %Vel perturbation m/s

% %--Point perturbation
dn=1; %10; %total width (grid points) of the perturbation (assumes an even number)
dn=(dn/2 -1);
   %z_final1(nzo-dn:nzo+dn,nxo-dn:nxo+dn)=c;   
z_final1(nzo,nxo)=c;   

%}
         
%%% Writing the model
        % %{   
       %file=['Vmodel_sigma20_pert_farfield_r5.xyz'];  %['test.xyz'];
       file=['VmodelElastic_Large_sigma20.xyz'];
      %    file=['VmodelElastic_Large_sigma5.xyz'];
               f=fopen(file,'w');
            [nygll,nxgll]=size(z_final); %size(zgll);
                for ix=1:nxgll
                    for iy=1:nygll
                        if iy>nygll-nygll1
                        vp=z_final1(iy,ix);  
                        vs=vp/sqrt(3);
                        elseif iy>(nygll-nygll2) & iy<(nygll-nygll1)
                        vp=z_final2(iy,ix);
                        vs=vp/sqrt(3);
                        else
                        vp=z_final1(iy,ix);  
                        vs=vp/sqrt(3);   
            
                        end
                    fprintf(f,'%10.2f %10.2f %10.2f %10.2f %10.2f\n',...
                   Xgll(iy,ix),Ygll(iy,ix),rho,vp,vs); %ORIGINAL 
                   % Xgll(iy,ix),Ygll(iy,ix),rho,vp,0); %ACOUSTIC vs=0
                    end
                end
            fclose(f);
            
            
            
            % CONTROL Vs>0 always (EL media)
            vs_min(count)=min(min([z_final1; z_final2] ));   % /sqrt(3)         
           
            
              %}

%% ----------Figures--------    
%{
      zfinal11=z_final1(1:nygll1,:);
zfinal22=z_final2(nygll1+1:nygll2,:);
zfinal33=z_final1(nygll2+1:end,:);
zfinal_r=[zfinal11 ;zfinal22; zfinal33];%% Dessins / Drawings



clims_fin=[vref*(1-sigma*3) vref*(1+sigma*3)]; %colormap(gray);
clims=[-sigma*3 sigma*3];

%Plot using Matlab Convention.
%In Specfem the model will be read Upside down
figure(20);
imagesc(1:xmax,1:ymax,zfinal_r,clims_fin); 
xlim([0 xmax]); ylim([0 ymax]);
colorbar;
title('zfinal(x,y)')
set(gca,'fontsize', 18);

figure(21);
imagesc(1:xmax,1:ymax,zfinal_r);  
xlim([0 xmax]); ylim([0 ymax]);
colorbar;


% %--Control--
vmax=max(max(zfinal_r));  %Vmax ondasP
vmin=min(min(zfinal_r));  %Vmin ondasP  vref*(1-sigma*4)
promax=abs(vref-vmax)*100/vref;
promin=abs(vref-vmin)*100/vref;
Vsmin=vmin/sqrt(3) % vmin S-waves > 0 if you want elastic medium!!

Vsmax=vmax/sqrt(3);  % I must satisfy Vsmax<vpmin  (Review this)
vmin-Vsmax  %vpmin-Vsmax>0
%{
std=std(zfinal_r(:));  %Standard deviation of Vp final  
mean(mean(zfinal_r));
pro_std=abs(std)*100/vref;  %-> Sigma (std deviation of the velocity in percentage!!)

figure(22);
histogram(zfinal_r(:))
%%%%%%%%histogram(z_final1(:))
title('Histogram of the P-VelModel')
xlabel('Vp (m/s)')
ylabel('Frequency')
set(gca,'fontsize', 18);

%}


%}
%Go back to the main directory
    %cd (['/Users/alejandro/Kernel_ElasticMedia_Project_Specfem/Matlab_Codes_Elastic/VARIOUS_Models2']); 
    cd ../..

        end
        
        
        %To control. All should be >0
vs_min=vs_min/sqrt(3);
%save vs_min  
        
toc


%{
% STATIONS  LARGE MEDIUM:
% Use only r: 1, 2, 5, 7, 15, 16, 17, 18 , 20, 21  of the original Medium in THIS  order
%The locations used in the paper so far are only r5, r18 and r20 (original medium)

%Coordinates of selected locations (Original Medium coordinates)
x_rp_o=[8400.24 8400.24 8400.24 8400.24 8400.24  8400.24  6533.52 5133.48 2177.84 933.36];
z_rp_o=[1594.49 2786.65 6300.18 8647.01 11005.87 16009.25 6300.18 6300.18 6300.18 6300.18];

%Coordinates in Large medium
xrcv=16800.48;  zrcv=16800.48-4705.69;  % My main anomaly is in the center of medium 
xrp_new=xrcv-(x_rp_o(1)-x_rp_o);
zrp_new=zrcv+abs(z_rp_o(1)-z_rp_o);
%}
