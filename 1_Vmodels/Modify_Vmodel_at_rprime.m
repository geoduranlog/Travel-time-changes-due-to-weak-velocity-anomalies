% Modify Vmodel at rprime (Pert Point) only
clear all; 
close all;


pert_percent=80; %7; %1.54; %100; %25;    %  vel anomaly of 7% at rprime
n=4;  %4 for modify P-waves   5 for S-waves only

%--- Set perturbation location----
% r'=r5 location (x,z)=(8400.24 , 6300.18)      z=0 is at the bottom

M=20 %1;
cd (['/Volumes/LaCie/Alejandro/KernelsComparison_paper_HARDDISK/Normalization_factor/OldSpecfem/Hetero_Elastic/Models/M',num2str(M,'%01.0f')]);


%Anomaly coordinates (m)
xo=8400.24; zo=6300.18;  %r5   %   8647.01;  %6300.18;  %3953.35;
%xo=5133.48; zo=6300.18; %r18
%xo=2177.84; zo=6300.18; %r20

%Read ASCII-delimited file of numeric data into matrix
Vmodel = dlmread('VmodelElastic_sigma20.xyz');

%Line of that point in the text file
line=find(   round(Vmodel(:,1)==xo) &  round(Vmodel(:,2)==zo)  );


%%  Create Pert Vel Models
for M=1:20% 20 %40 %nmodels
    %M=1
    
    %%Go to the model directory 
    cd (['/Volumes/LaCie/Alejandro/KernelsComparison_paper_HARDDISK/Normalization_factor/OldSpecfem/Hetero_Elastic/Models/M',num2str(M,'%01.0f')]);
    
    

%Read ASCII-delimited file of numeric data into matrix
Vmodel = dlmread('VmodelElastic_sigma20.xyz');


v_rprime=Vmodel(line, n);  %P or S-vel at r'
vpert=v_rprime+v_rprime*pert_percent/100;  
%vpert=v_rprime+100;  %Weak Anomaly in Oberman 2013 vpert=v+100m/s

%Modify Vel Model at r' location - One waves type only
%Vmodel(line, n)=vpert;  

%In case of Both P & S.  OJO! Choose v_rprime for P-waves (n=4) and Then:
Vmodel(line, 4)=vpert;   Vmodel(line, 5)=vpert/sqrt(3); 

  
%format shortE   %Standard is 'format short'  How show values


%Write New model    
%cd (['/Users/Alejandro/KernelsComparison_paper/Normalization_factor/OldSpecfem/Hetero_Elastic/Models/M',num2str(M,'%01.0f')]);
cd (['/Volumes/LaCie/Alejandro/KernelsComparison_paper_HARDDISK/Normalization_factor/OldSpecfem/Hetero_Elastic/Models/M',num2str(M,'%01.0f')]);
    

            %file=['VmodelEl_sigma20_Vs_pert7pro_r5.xyz'];  
                  % file=['VmodelEl_sigma20_Vp_pert7pro_r2.xyz'];  
                  % file=['VmodelEl_sigma20_pert25pro_r5.xyz'];  
         %   file=['VmodelEl_sigma20_Vs__pert100pro_r20.xyz']; 
       %file=['VmodelEl_sigma20_pert100pro_r20.xyz'];  %Weak was actually        %1.75% . I wantes 1.54%
       file=['VmodelEl_sigma20_pert80pro_r5.xyz'];   
       
            f=fopen(file,'w');
            fprintf(f,'%10.2f %10.2f %10.2f %10.2f %10.2f\n',Vmodel')
                      fclose(f);

                                          
                      end