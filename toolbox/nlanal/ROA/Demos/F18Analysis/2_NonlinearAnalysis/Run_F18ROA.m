%
% This file runs the ROA V-S Iteration and MC Sim for F/A-18 
%
% Abhijit 03/05/2010

%---------------------------------------------------------

LB = 1;  UB = 0; 
d2r = pi/180; r2d = 1/ d2r; 
% --------------------------------------------------------------------
% Run Revised

load 2_Mar14_PolyF18RevisedModel_Phi_35
f = cleanpoly(xcldotR,10^-06);

if LB
    V0 = strcat('1_Mar14_QuadROAVS_F18Revised_',num2str(xtrim(7)*r2d));
    ROA2_F18_VSIteration; 
elseif UB
    V0 = strcat('1_Mar14_ROAMC_F18Revised_',num2str(xtrim(7)*r2d));
    Search4UnstableTraj;
end


% --------------------------------------------------------------------
% Run Baseline 

load 2_Mar14_PolyF18BaseLineModel_Phi_35
f = cleanpoly(xcldotB,10^-06);  
 

if LB
    V0 = strcat('1_Mar14_QuadROAVS_F18Baseline_',num2str(xtrim(7)*r2d));
    ROA2_F18_VSIteration;
elseif UB
    V0 = strcat('1_Mar14_ROAMC_F18BaseLine_',num2str(xtrim(7)*r2d));
    Search4UnstableTraj;
end


