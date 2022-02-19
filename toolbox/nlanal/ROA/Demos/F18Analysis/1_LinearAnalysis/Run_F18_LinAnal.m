%==========================================================================
% Generate Linear Plants
if steady
    % --------  Coordinated Case
    Bank0  = load('F18Model_Coord_Phi_0', 'Pss','CssB','CssR') ;   
    Bank10 = load('F18Model_Coord_Phi_10','Pss'); 
    Bank25 = load('F18Model_Coord_Phi_25','Pss'); 
    Bank35 = load('F18Model_Coord_Phi_35','Pss'); 
    
else
    % --------  UnCoordinated Case
    Bank0  = load('F18Model_UnCoord_Phi_0','Pss','CssB','CssR') ;   
    Bank10 = load('F18Model_UnCoord_Phi_10','Pss'); 
    Bank25 = load('F18Model_UnCoord_Phi_25','Pss'); 
    Bank35 = load('F18Model_UnCoord_Phi_35','Pss'); 

end

Pss_Bank0 =  Bank0.Pss;   
Pss_Bank10 = Bank10.Pss;  
Pss_Bank25 = Bank25.Pss;   
Pss_Bank35 = Bank35.Pss;   

Css_Baseline = Bank0.CssB;         
Css_Revised =  Bank0.CssR;

%-----------------------------------------------------------------------
% Stack The Stable Closed Loop LTI Models
              
Pss_Plants(:,:,1:4) = stack(1, Pss_Bank0, Pss_Bank10, Pss_Bank25,Pss_Bank35); 
             
            
%==========================================================================
%                       Linear Analysis
%==========================================================================



%==========================================================================
% Loopmargin Analysis
% 
%==========================================================================

smiB = cell(4,1);  smiR = cell(4,1); 
dmiB = cell(4,1);  dmiR = cell(4,1); 
mmiB = cell(4,1);  mmiR = cell(4,1); 
dgb  = cell(4,1);  dpb = cell(4,1); 
dgr  = cell(4,1);  dpr = cell(4,1); 

for i1 = 1:4 
    
    %--------------------------------------------------------------------------
    % Baseline Control Law
    [smiB{i1,1},dmiB{i1,1},mmiB{i1,1},smoB,dmoB,mmoB,mmioB]= ...
          loopmargin(Pss_Plants(:,:,i1,1),Css_Baseline);


    %--------------------------------------------------------------------------
    % Revised Control Law
    [smiR{i1,1},dmiR{i1,1},mmiR{i1,1},smoR,dmoR,mmoR,mmioR]=...
                 loopmargin(Pss_Plants(:,:,i1,1),Css_Revised);

    % Extract the DMPLOT results to plot Disk Margin  
    [dgb{i1,1},dpb{i1,1}] = dmplot(max(mmiB{i1,1}.GainMargin));
    [dgr{i1,1},dpr{i1,1}] = dmplot(max(mmiR{i1,1}.GainMargin)); 

end
%--------------------------------------------------------------------------
% Plot Results

for i1 = 1:4
    figure(i1)
    ph1 = plot(dgb{i1,:},dpb{i1,:},'--r');
    hold on 
    ph2 = plot(dgr{i1,:},dpr{i1,:},'-b'); 
    legend('Baseline','Revised','Location','Best')
    xl=xlabel('Gain Margin (dB)'); 
    yl=ylabel('Phase Margin (deg)');
    set(ph1,'LineWidth',2); set(ph2,'LineWidth',2);
    set([xl,yl,gca],'FontSize',16);
    grid on 
    hold off
end



disp('End Loopmargin Analysis')

return
%==========================================================================


%==========================================================================
% Input Multiplicative Uncertainty Analysis
% (a)  Diagonal Input Multiplicative  Uncertainty 
% (b)  Full Block Input Multiplicative  Uncertainty
%==========================================================================


%==========================================================================
% Diagonal Input Multiplicative  Uncertainty 

d1=ultidyn('d1');
d2=ultidyn('d2');
d3=ultidyn('d3');
DeltaDiag= diag([d1,d2,d3]);

%--------------------------------------------------------------------------
% Frequency shaping of error dynamics

% Frequency Shaping of weighting function.  We have chosen it to be 1. 
W_in = eye(3); 


for i1 = 1:4
    
   
    %--------------------------------------------------------------------------
    % Building up Baseline Uncertain Model : Diagonal

    PuncB = Pss_Plants(:,:,i1,1)*(eye(3)+W_in*DeltaDiag);

    % Uncertain Baseline Closed loop
    CLuncB = feedback(PuncB,Css_Baseline);

    om = logspace(-2,2,100); 
    CLuncB = frd(CLuncB,om); 

    %--------------------------------------------------------------------------
    % Robust Stability Analysis for Baseline Model : Diagonal
    [STABMARGB,DESTABUNCB,REPORTB,INFOB] = robuststab( CLuncB );
    BaselineMUDiag(:,i1) = INFOB.MussvBnds(1); 

    %--------------------------------------------------------------------------
    % Building up Revised Uncertain Model  : Diagonal  

    PuncR = Pss_Plants(:,:,i1,1)*(eye(3)+W_in*DeltaDiag);

    % Uncertain Revised Closed loop
    CLuncR = feedback(PuncR,Css_Revised);

    %om = logspace(-2,2,100); 
    CLuncR = frd(CLuncR,om); 

    %--------------------------------------------------------------------------
    % Robust Stability Analysis for Revised Model : Diagonal
    [STABMARGR,DESTABUNCR,REPORTR,INFOR] = robuststab( CLuncR );
    RevisedMUDiag(:,i1) = INFOR.MussvBnds(1);  

end

%--------------------------------------------------------------------------
% Plot Robustness Results
    
figure(2)    
ph1 = semilogx(BaselineMUDiag(:,1:4),'--r'); 
legend(ph1, 'Baseline')
hold on
ph2 = semilogx(RevisedMUDiag(:,1:4),'-b');
xl=xlabel('frequency (rad/s)'); 
yl=ylabel('\mu');
set(ph1,'LineWidth',2); set(ph2,'LineWidth',2);
set([xl,yl,gca],'FontSize',16);
grid on 
legend([ph1(1) ph2(1)],'Baseline','Revised') 


%==========================================================================
% Full Block Input Multiplicative Uncertainty 
        
DeltaFull = ultidyn('Delta',[3 3],'Bound',1);


for i1 = 1:4
    
    %--------------------------------------------------------------------------
    % Building up Baseline Uncertain Model : Full Block Uncertainty

    PuncB = Pss_Plants(:,:,i1,1)*(eye(3)+W_in*DeltaFull);

    % Uncertain Baseline Closed loop
    CLuncB = feedback(PuncB,Css_Baseline);

    om = logspace(-2,2,100); 
    CLuncB = frd(CLuncB,om); 

    %--------------------------------------------------------------------------
    % Robust Stability Analysis for Baseline Model: Full Block Uncertainty
    [STABMARGB,DESTABUNCB,REPORTB,INFOB] = robuststab( CLuncB );
    BaselineMUFull(:,i1) = INFOB.MussvBnds(1); 



    %--------------------------------------------------------------------------
    % Building up Revised Uncertain Model  : Full Block Uncertainty  

    PuncR = Pss_Plants(:,:,i1,1)*(eye(3)+W_in*DeltaFull);

    % Uncertain Revised Closed loop
    CLuncR = feedback(PuncR,Css_Revised);

    %om = logspace(-2,2,100); 
    CLuncR = frd(CLuncR,om); 
    %--------------------------------------------------------------------------
    % Robust Stability Analysis for Revised Model: Full Block Uncertainty 
    [STABMARGR,DESTABUNCR,REPORTR,INFOR] = robuststab( CLuncR );
    RevisedMUFull(:,i1) = INFOR.MussvBnds(1);  

end
%--------------------------------------------------------------------------
% Plot Robustness Results
    
figure(3)
ph1 = semilogx(BaselineMUFull(:,1:4),'--r');
hold on
ph2 = semilogx(RevisedMUFull(:,1:4),'-b');
xl=xlabel('frequency (rad/s)'); 
yl=ylabel('\mu');
set(ph1,'LineWidth',2); set(ph2,'LineWidth',2);
set([xl,yl,gca],'FontSize',16);
grid on 
legend([ph1(1) ph2(1)],'Baseline','Revised') 


disp('End Input Multiplicative Uncertainty Analysis')

%==========================================================================




%==========================================================================
% Parametric Uncertainty in Aerodynamic Coefficient Analysis
%==========================================================================

    
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
var_unc = 10; 

for i1 = 1:4
    
    Pss = Pss_Plants(:,:,i1,1); 
    

    %==========================================================================
    % Analysis on  Baseline Model  

    % Declare Uncertain Variable : Real Uncertainty
    d1 = ureal('d1',Pss.A(1,1),'pe',var_unc); 
    d2 = ureal('d2',Pss.A(3,1),'pe',var_unc);
    d3 = ureal('d3',Pss.A(5,1),'pe',var_unc); 
    d4 = ureal('d4',Pss.A(5,5),'pe',var_unc);
    d5 = ureal('d5',Pss.A(3,3),'pe',var_unc); 
    d6 = ureal('d6',Pss.A(2,4),'pe',var_unc); 
    d7 = ureal('d7',Pss.A(4,2),'pe',var_unc); 
    d8 = ureal('d8',Pss.A(4,4),'pe',var_unc);


    %--------------------------------------------------------------------------
    % Uncertain A Matrix
    Aunc = [ d1 Pss.A(1,2:end); ...
             Pss.A(2,1:3) d6 Pss.A(2,5:6);...
            d2 Pss.A(3,2) d5 Pss.A(3,4:end); ...
            Pss.A(4,1) d7 Pss.A(4,3) d8 Pss.A(4,5:6);...
            d3 Pss.A(5,2:4) d4 Pss.A(5,6); ...
            Pss.A(6,:) ]; 

    %--------------------------------------------------------------------------
    % Uncertain Open Loop Baseline Model 
    PuncB  = ss(Aunc,Pss.B,Pss.C,Pss.D);

    %--------------------------------------------------------------------------
    % Uncertain Closed loop
    CLuncB = feedback(PuncB,Css_Baseline);


    %--------------------------------------------------------------------------
    % Robust Stability / Mu Analysis for Baseline Model

    CLuncB = complexify(CLuncB,0.05);

    om = logspace(-2,2,100); 
    CLuncB = frd(CLuncB,om); 

    opt2=robopt('Mussv','m9');
    [SMB , DeUncB, repB, infoB] = robuststab(CLuncB,opt2);
    BaselinePara(:,i1) = infoB.MussvBnds(1);



    %==========================================================================
    % Analysis on Revised Model  

    %--------------------------------------------------------------------------
    % Declare Uncertain Variable : Real Uncertainty
    d1 = ureal('d1',Pss.A(1,1),'pe',var_unc); 
    d2 = ureal('d2',Pss.A(3,1),'pe',var_unc);
    d3 = ureal('d3',Pss.A(5,1),'pe',var_unc); 
    d4 = ureal('d4',Pss.A(5,5),'pe',var_unc);
    d5 = ureal('d5',Pss.A(3,3),'pe',var_unc); 
    d6 = ureal('d6',Pss.A(2,4),'pe',var_unc); 
    d7 = ureal('d7',Pss.A(4,2),'pe',var_unc); 
    d8 = ureal('d8',Pss.A(4,4),'pe',var_unc);



    %--------------------------------------------------------------------------
    % Uncertain A Matrix
    Aunc = [ d1 Pss.A(1,2:end); ...
             Pss.A(2,1:3) d6 Pss.A(2,5:6);...
             d2 Pss.A(3,2) d5 Pss.A(3,4:end); ...
             Pss.A(4,1) d7 Pss.A(4,3) d8 Pss.A(4,5:6);...
             d3 Pss.A(5,2:4) d4 Pss.A(5,6); ...
             Pss.A(6,:) ]; 


    %--------------------------------------------------------------------------
    % Uncertain Open Loop Revised Model 
    PuncM  = ss(Aunc,Pss.B,Pss.C,Pss.D);

    %--------------------------------------------------------------------------
    % Uncertain  Revised Closed loop
    CLuncM = feedback(PuncM,Css_Revised);


    %--------------------------------------------------------------------------
    % Robust Stability / Mu Analysis

    CLuncM = complexify(CLuncM,0.05);

    %om = logspace(-2,3,100); 
    CLuncM = frd(CLuncM,om); 

    opt2=robopt('Mussv','m9');
    [SMM , DeUncM, repM, infoM] = robuststab(CLuncM,opt2);
    RevisedPara(:,i1) = infoM.MussvBnds(1);
end    
%--------------------------------------------------------------------------
% Plot Results

figure(4)
ph1 = semilogx(BaselinePara(:,1:4),'-.r');
hold on
ph2 = semilogx(RevisedPara(:,1:4),'-b');
xl=xlabel('frequency (rad/s)'); 
yl=ylabel('\mu');
set(ph1,'LineWidth',2); set(ph2,'LineWidth',2);
set([xl,yl,gca],'FontSize',16);
grid on 
legend([ph1(1) ph2(1)],'Baseline','Revised') 
% End Parametric Uncertainty in Aerodynamic Coefficient Analysis
%==========================================================================
