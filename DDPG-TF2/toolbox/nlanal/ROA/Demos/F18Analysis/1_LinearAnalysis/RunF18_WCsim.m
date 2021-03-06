%
% This file perfroms a worst-case analysis of the F/A-18 linear plants. 

close all
steady = 0; 


%==========================================================================
% Generate Nominal Linear Plants
if steady
    % --------  Coordinated Case
    Bank0  = load('F18Model_Coord_Phi_0', 'Pss','CssB','CssR') ;   
    Bank10 = load('F18Model_Coord_Phi_10','Pss'); 
    Bank25 = load('F18Model_Coord_Phi_25','Pss'); 
    Bank35 = load('F18Model_Coord_Phi_35','Pss','CssB','CssR'); 
    
else
    % --------  UnCoordinated Case
    Bank0  = load('F18Model_UnCoord_Phi_0','Pss','CssB','CssR') ;   
    Bank10 = load('F18Model_UnCoord_Phi_10','Pss'); 
    Bank25 = load('F18Model_UnCoord_Phi_25','Pss'); 
    Bank35 = load('F18Model_UnCoord_Phi_35','Pss','CssB','CssR'); 

end

Pss_Bank0 =  Bank0.Pss;   
Pss_Bank10 = Bank10.Pss;  
Pss_Bank25 = Bank25.Pss;   
Pss_Bank35 = Bank35.Pss;   

Css_Baseline = Bank0.CssB;         
Css_Revised =  Bank0.CssR;

%----------------------------------------
% Stack The Stable Closed Loop LTI Models
              
Pss_Plants(:,:,1:4) = stack(1, Pss_Bank0, Pss_Bank10, Pss_Bank25,Pss_Bank35); 


%==========================================================================
% Construct Uncertain Plant
%==========================================================================
%
% 1.  parametric uncertainty 
% 2.  Diagonal IMU unmodeled dynamics
% 3.  Full Block IMU


%----------------------------------------------------------
% Parametric Uncertainty in Aerodynamic Coefficient Analysis
%----------------------------------------------------------
    
% Lateral Direction Important Terms

%  Uncertain Entries    Remarks 
%  in A Matrix
%  (1,1)                Sidefore due to sideslip      
%  (3,1)                Rolling moment due to sideslip
%  (5,1)                Yawing moment due to sideslip
%  (5,5)                Yaw Damping
%  (3,3)                Roll Damping

% Longitudinal Direction Important Terms

%  Uncertain Entries    Remarks   
%  in A Matrix
%  (2,4)                Normal fore due to Pitch Rate      
%  (4,2)                Pitch Stiffness
%  (4,4)                Pitch Damping

% Uncertainty Variation in the real parameters
var_unc  = 10; 


%----------------------------------------------------------
% IMU Uncertainty in 
%----------------------------------------------------------

DeltaFull = ultidyn('DeltaFull',[3 3],'Bound',1);

d11     = ultidyn('d11');
d22     = ultidyn('d22');
d33     = ultidyn('d33');
Delta   = diag([d11,d22,d33]);


% Initialize Model Storage 
Punc      = cell(4,1); 
PuncPara  = cell(4,1); 
PuncDFull = cell(4,1); 
PuncDdiag = cell(4,1); 
% CLuncB    = cell(4,1); 
% CLuncR    = cell(4,1); 



for i1 = 1:4
    
    Pss = Pss_Plants(:,:,i1,1); 
    

    %===========================================================
    % Analysis on Model  

    % Declare Uncertain Variable : Real Uncertainty
    d1 = ureal('d1',Pss.A(1,1),'pe',var_unc); 
    d2 = ureal('d2',Pss.A(3,1),'pe',var_unc);
    d3 = ureal('d3',Pss.A(5,1),'pe',var_unc); 
    d4 = ureal('d4',Pss.A(5,5),'pe',var_unc);
    d5 = ureal('d5',Pss.A(3,3),'pe',var_unc); 
    d6 = ureal('d6',Pss.A(2,4),'pe',var_unc); 
    d7 = ureal('d7',Pss.A(4,2),'pe',var_unc); 
    d8 = ureal('d8',Pss.A(4,4),'pe',var_unc);


    %------------------------------------------------------------
    % Uncertain A Matrix
    Aunc = [ d1 Pss.A(1,2:end); ...
            Pss.A(2,1:3) d6 Pss.A(2,5:6);...
            d2 Pss.A(3,2) d5 Pss.A(3,4:end); ...
            Pss.A(4,1) d7 Pss.A(4,3) d8 Pss.A(4,5:6);...
            d3 Pss.A(5,2:4) d4 Pss.A(5,6); ...
            Pss.A(6,:) ]; 

    %--------------------------------------------------------------
    % Uncertain Open Loop Model 
    P1 = ss(Aunc,Pss.B,Pss.C,Pss.D);
    scl = 1000; 
    PuncPara{i1,1}    = P1;
    PuncDFull{i1,1}   = (P1.nom/scl)*(eye(3)+ DeltaFull);
    PuncDdiag{i1,1}   = (P1.nom/scl)*(eye(3)+ Delta);
%     CLuncB{i1,1} = feedback(Punc{i1,1},Css_Baseline); 
%     CLuncR{i1,1} = feedback(Punc{i1,1},Css_Revised); 
  
end



%==========================================================================
% Run Worst Case Analysis
%==========================================================================

% ----------------------------------------------------------
% WC Analysis on Paramteric uncertainty
% ---------------------------------------------------------




for j1 = 3
    
        if j1 == 1
             P1 = PuncPara{1,1};
             P2 = PuncPara{2,1};
             P3 = PuncPara{3,1};
             P4 = PuncPara{4,1};
             scl = 1; 
        elseif j1 == 2
             P1 = PuncDFull{1,1};
             P2 = PuncDFull{2,1};
             P3 = PuncDFull{3,1};
             P4 = PuncDFull{4,1};
        elseif j1 == 3
             P1 =  PuncDdiag{1,1};
             P2 =  PuncDdiag{2,1};
             P3 =  PuncDdiag{3,1};
             P4 =  PuncDdiag{4,1};
        end

        % -------- Nominal Analysis
        LTnom_B1 = loopsens(P1,Css_Baseline); 
        LTnom_R1 = loopsens(P1,Css_Revised); 
        
        LTnom_B2 = loopsens(P2,Css_Baseline); 
        LTnom_R2 = loopsens(P2,Css_Revised);
        
        LTnom_B3 = loopsens(P3,Css_Baseline); 
        LTnom_R3 = loopsens(P3,Css_Revised);
        
        LTnom_B4 = loopsens(P4,Css_Baseline); 
        LTnom_R4 = loopsens(P4,Css_Revised);
        
        


        % ------- Pointwise Frequency
        opt = wcgopt('FreqPtWise',1);
        %opt = wcgopt; 

        % ------- Frequency gridding
        omega = logspace(-2,2,50); 
        LTnom_BSiBetaFrd1 = frd(LTnom_B1.PSi(5,1:2),omega); 
        LTnom_RSiBetaFrd1 = frd(LTnom_R1.PSi(5,1:2),omega); 

        LTnom_BSiBetaFrd2 = frd(LTnom_B2.PSi(5,1:2),omega); 
        LTnom_RSiBetaFrd2 = frd(LTnom_R2.PSi(5,1:2),omega);
        
        LTnom_BSiBetaFrd3 = frd(LTnom_B3.PSi(5,1:2),omega); 
        LTnom_RSiBetaFrd3 = frd(LTnom_R3.PSi(5,1:2),omega);
        
        LTnom_BSiBetaFrd4 = frd(LTnom_B4.PSi(5,1:2),omega); 
        LTnom_RSiBetaFrd4 = frd(LTnom_R4.PSi(5,1:2),omega);
        % 
        % ------- Worst Case Gain 
        [maxgainBBeta1,maxgainuncBBeta1,infoBBeta1] = wcgain(LTnom_BSiBetaFrd1,opt);
        [maxgainRBeta1,maxgainuncRBeta1,infoRBeta1] = wcgain(LTnom_RSiBetaFrd1,opt);
        
        [maxgainBBeta2,maxgainuncBBeta2,infoBBeta2] = wcgain(LTnom_BSiBetaFrd2,opt);
        [maxgainRBeta2,maxgainuncRBeta2,infoRBeta2] = wcgain(LTnom_RSiBetaFrd2,opt);
        
        [maxgainBBeta3,maxgainuncBBeta3,infoBBeta3] = wcgain(LTnom_BSiBetaFrd3,opt);
        [maxgainRBeta3,maxgainuncRBeta3,infoRBeta3] = wcgain(LTnom_RSiBetaFrd3,opt);
        
        [maxgainBBeta4,maxgainuncBBeta4,infoBBeta4] = wcgain(LTnom_BSiBetaFrd4,opt);
        [maxgainRBeta4,maxgainuncRBeta4,infoRBeta4] = wcgain(LTnom_RSiBetaFrd4,opt);
        %---------- Plot Results
        
      
        
        figure 
        set(gca,'FontSize',20)
        
        semilogx(scl*maxgainBBeta1.UpperBound,'--r',scl*maxgainRBeta1.UpperBound,'-b');
         hold on
        semilogx(scl*maxgainBBeta2.UpperBound,'--r',scl*maxgainRBeta2.UpperBound,'-b');
        semilogx(scl*maxgainBBeta3.UpperBound,'--r',scl*maxgainRBeta3.UpperBound,'-b');
        semilogx(scl*maxgainBBeta4.UpperBound,'--r',scl*maxgainRBeta4.UpperBound,'-b');
        xlabel('Frequency (rad/sec)','FontSize',20)
        ylabel('Worst Closed-loop Gain','FontSize',20);
        legend('Baseline','Revised','Location','NorthEast');
        h = findobj(gcf,'type','line'); set(h,'linewidth',2);
        grid on 
        hold off
        
%         figure
%         hold off
%         semilogx(fnorm(LTnom_BSiBetaFrd.nom),'-b')
%         hold on 
%         semilogx(fnorm(LTnom_RSiBetaFrd.nom),'-r')
%         legend('Baseline','Revised')
%         xlabel('Frequency (rad/sec)')
%         ylabel('Nominal Closed-loop Gain of d to e');
%         legend('Baseline','Revised','Location','NorthWest');
%         hold off
end

