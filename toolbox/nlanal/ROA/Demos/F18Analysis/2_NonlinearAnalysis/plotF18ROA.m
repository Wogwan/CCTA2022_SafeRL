% Plot ROA Results:


%===================================================================

d2r = pi/180;
r2d = 1/d2r;

%===================================================================
% Shape Matrix (Un-Normalized)
% x     =     [beta  ;  alpha;   p;        q;     r;       phi;      xc];
% Dmax = diag([5*d2r    25*d2r   25*d2r    25*d2r 5*d2r   45*d2r   20*d2r]);  

%==========================================================================
% Create Shape Function
d2r = pi/180;
Dmax = diag([10*d2r 25*d2r 35*d2r 30*d2r 15*d2r 25*d2r 20*d2r]);  
N = inv(Dmax^2);   %  Scale by inverse of max state values
%N = N/max(N(:));   %  Normalize by the largest entry

pvar beta alpha p q r phi xc
x = [ beta; alpha; p; q; r ; phi; xc]; 
p = x'*N*x; 

%===================================================================
%  Alpha - Beta 

Pba = diag([10*d2r 25*d2r]);     %  Alpha - Beta 
N1 = inv(Pba^2);                %  Scale by inverse of max state values
Pba = N1;%/max(N1(:));            %  Normalize by the largest entry
Prt1 = sqrtm(Pba);


%===================================================================
% p-r
Ppr = diag([35*d2r  15*d2r]);     % p-r
N2 = inv(Ppr^2);                 %  Scale by inverse of max state values
Ppr = N2;%/max(N2(:));             %  Normalize by the largest entry
Prt2 = sqrtm(Ppr);
 

%===================================================================
% 35 DEGREE BANK ANGLE 

ftrim = load('2_Mar14_PolyF18RevisedModel_Phi_35','xtrim');

alphatrim = ftrim.xtrim(3);      betatrim = ftrim.xtrim(2);
ptrim     = ftrim.xtrim(4);      rtrim    = ftrim.xtrim(6); 


 BClawVs = load('1_Mar14_ROAVS_F18Baseline_35');
 RClawVs = load('1_Mar14_ROAVS_F18Revised_35'); 
% p(x) = beta Value:
% beta_ROA =  min(info.BetaBnds); 


% --------------------Final ROA Vs Result
BbetaVs =  BClawVs.info.iteration(end).beta;
RbetaVs =  RClawVs.info.iteration(end).beta;

% --------------------Linearized ROA Vs Result
BbetaVslin = BClawVs.info.iteration(1).beta;
RbetaVslin = RClawVs.info.iteration(1).beta;

% -------------------- Quad ROA Vs Result
 BClawVsQuad = load('1_Mar14_QuadROAVS_F18Baseline_35');
 RClawVsQuad = load('1_Mar14_QuadROAVS_F18Revised_35'); 
 BbetaVsQuad = BClawVsQuad.info.iteration(end).beta;
 RbetaVsQuad = RClawVsQuad.info.iteration(end).beta;


% --------------------Monte Carlo Result
 BClawMC = load('1_Mar14_ROAMC_F18Baseline_35');
 RClawMC = load('1_Mar14_ROAMC_F18Revised_35'); 

% simply scaling the IC found from simulation has been able to reduce the
% MC beta value

% x01 = BClawMC.xun(:,i); x01 = x01 / 1.0386; 
% x01 = min(RClawMC.betaun);  x01 = x01 / 1.0206 

BbetaMC =  min(BClawMC.betaun)* 0.985;   % Smaller unstable traj was found by hand
RbetaMC =  min(RClawMC.betaun);  

%===================================================================
% Gridding Data
Nth = 100;
th = linspace(0,2*pi,Nth);


%===================================================================
% Revised Simulation 
% level Set calculation

yy_ROA = [cos(th); sin(th)];

abtrim = [betatrim*ones(1,Nth); alphatrim*ones(1,Nth)];
prtrim = [ptrim*ones(1,Nth); rtrim*ones(1,Nth)];
%--------------------------------------------------
%  alpha - beta

xxba_ROA_V4_B = sqrt(BbetaVs)*inv(Prt1)*yy_ROA + abtrim  ; 
xxba_ROA_V4_R = sqrt(RbetaVs)*inv(Prt1)*yy_ROA + abtrim; 

xxba_ROA_V4_Blin = sqrt(BbetaVslin)*inv(Prt1)*yy_ROA + abtrim  ; 
xxba_ROA_V4_Rlin = sqrt(RbetaVslin)*inv(Prt1)*yy_ROA + abtrim; 

% xxba_ROA_V2_B = sqrt(BbetaVsQuad)*inv(Prt1)*yy_ROA + abtrim  ; 
% xxba_ROA_V2_R = sqrt(RbetaVsQuad)*inv(Prt1)*yy_ROA + abtrim; 

xxba_ROA_MC_B = sqrt(BbetaMC)*inv(Prt1)*yy_ROA + abtrim; 
xxba_ROA_MC_R = sqrt(RbetaMC)*inv(Prt1)*yy_ROA + abtrim; 

%--------------------------------------------------
%  p - r

xxpr_ROA_V4_B = sqrt(BbetaVs)*inv(Prt2)*yy_ROA + prtrim; 
xxpr_ROA_V4_R = sqrt(RbetaVs)*inv(Prt2)*yy_ROA + prtrim ; 

xxpr_ROA_V4_Blin = sqrt(BbetaVslin)*inv(Prt2)*yy_ROA + prtrim; 
xxpr_ROA_V4_Rlin = sqrt(RbetaVslin)*inv(Prt2)*yy_ROA + prtrim ; 

% xxpr_ROA_V2_B = sqrt(BbetaVsQuad)*inv(Prt2)*yy_ROA + prtrim; 
% xxpr_ROA_V2_R = sqrt(RbetaVsQuad)*inv(Prt2)*yy_ROA + prtrim ; 

xxpr_ROA_MC_B = sqrt(BbetaMC)*inv(Prt2)*yy_ROA + prtrim ; 
xxpr_ROA_MC_R = sqrt(RbetaMC)*inv(Prt2)*yy_ROA + prtrim ;

%===================================================================
% Plot Results
close all 


    
figure1 = figure;
% Create axes
axes('Parent',figure1,'FontSize',16) 
plot(xxba_ROA_V4_B(2,:)*r2d,xxba_ROA_V4_B(1,:)*r2d,'-b','LineWidth',2);
hold on
plot(xxba_ROA_MC_B(2,:)*r2d,xxba_ROA_MC_B(1,:)*r2d,'--r','LineWidth',1.5);
hold on
legend('Quartic Lyapunov','Monte Carlo Upper Bound')
plot(xxba_ROA_V4_R(2,:)*r2d,xxba_ROA_V4_R(1,:)*r2d,'-b','LineWidth',2);
hold on 
plot(xxba_ROA_MC_R(2,:)*r2d,xxba_ROA_MC_R(1,:)*r2d,'--r','LineWidth',1.5);
hold on 
plot(alphatrim*r2d,betatrim*r2d,'+','MarkerSize',12,'LineWidth',2)
grid on
ylabel('\beta (deg)','FontSize',16)
xlabel('\alpha (deg)','FontSize',16)

figure1 = figure;
% Create axes
axes('Parent',figure1,'FontSize',16)
plot(xxpr_ROA_V4_B(2,:)*r2d,xxpr_ROA_V4_B(1,:)*r2d,'-b','LineWidth',2);
hold on
plot(xxpr_ROA_MC_B(2,:)*r2d,xxpr_ROA_MC_B(1,:)*r2d,'--r','LineWidth',2);
legend('Quartic Lyapunov','Monte Carlo Upper Bound')
plot(xxpr_ROA_V4_R(2,:)*r2d,xxpr_ROA_V4_R(1,:)*r2d,'-b','LineWidth',2);
plot(xxpr_ROA_MC_R(2,:)*r2d,xxpr_ROA_MC_R(1,:)*r2d,'--r','LineWidth',2); 
plot(rtrim*r2d, ptrim*r2d,'+','MarkerSize',12,'LineWidth',2)
grid on
ylabel('p (deg/s)','FontSize',16)
xlabel('r (deg/s)','FontSize',16)


% figure1 = figure;
% % Create axes
% axes('Parent',figure1,'FontSize',16)
% % plot(xxpr_ROA_V2_B(2,:)*r2d,xxpr_ROA_V2_B(1,:)*r2d,'-b','LineWidth',2);
% % hold on
% % plot(xxpr_ROA_V2_R(2,:)*r2d,xxpr_ROA_V2_R(1,:)*r2d,'--r','LineWidth',2);
% hold on 
% plot(xxpr_ROA_V4_Blin(2,:)*r2d,xxpr_ROA_V4_Blin(1,:)*r2d,'-b','LineWidth',2);
% plot(xxpr_ROA_V4_Rlin(2,:)*r2d,xxpr_ROA_V4_Rlin(1,:)*r2d,'--r','LineWidth',2);
% plot(rtrim*r2d, ptrim*r2d,'+','MarkerSize',12,'LineWidth',2)
% legend('Baseline','Revised')
% grid on
% ylabel('p (deg/s)','FontSize',16)
% xlabel('r (deg/s)','FontSize',16)
% 
% figure1 = figure;
% % Create axes
% axes('Parent',figure1,'FontSize',16)
% % plot(xxba_ROA_V2_B(2,:)*r2d,xxba_ROA_V2_B(1,:)*r2d,'-b','LineWidth',2);
% % hold on
% % plot(xxba_ROA_V2_R(2,:)*r2d,xxba_ROA_V2_R(1,:)*r2d,'--r','LineWidth',2);
% hold on 
% plot(xxba_ROA_V4_Blin(2,:)*r2d,xxba_ROA_V4_Blin(1,:)*r2d,'-b','LineWidth',2);
% plot(xxba_ROA_V4_Rlin(2,:)*r2d,xxba_ROA_V4_Rlin(1,:)*r2d,'--r','LineWidth',2);
% plot(alphatrim*r2d,betatrim*r2d,'+','MarkerSize',12,'LineWidth',2)
% legend('Baseline','Revised')
% grid on
% ylabel('\beta (deg)','FontSize',16)
% xlabel('\alpha (deg)','FontSize',16)




return
