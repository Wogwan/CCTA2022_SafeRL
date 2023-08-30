%
% This file compares the Closed  Loop Models  
%

% ===================================================================
% ------- Full Nonlinear Model

% Cases for Trimming: 
% 
Level = 0; 
CoordTurn = 1;
UnCoordTurn = 0;

% --------- SET UP Full Nonlinear CLP Sim Setting 

displayon   = 1;          % Display Trim values
opts.poly   = 0; 

GenF18LinModel; 
GenF18LinCLPModel; 

% ------- Sim paramter setup
xtf = xtrim;  utf = utrim; 
C = Pss.C;  D = Pss.D; 
aytrim = double(subs(C(1,:)*x + D(1,:)*u,[x;u],[xtf(2:7);utf(1:3)]));
bdottrim = double(subs(C(end,:)*x + D(end,:)*u,[x;u],[xtf(2:7);utf(1:3)]));

InpPert = 0;   % Perturb the input channel 
dxic = [0; 5; 20; 25 ; 20; 20; 5; 0; 0]*d2r;  % initial cond on simulation 
dxic = [0; 10; 20; 25 ; 25; 20; 25; 0; 0]*d2r;
dx0   = xtf+dxic ;
Tfinal = 15; 

% ---- Baseline Simulation 
K = CssB; 
[tB, yBresp] = sim('f18nlclp',[0 Tfinal]); 


% ---- Revised Simulation 

    
K = CssR; 
[tR, yRresp] = sim('f18nlclp',[0 Tfinal]); 



% ===================================================================
% ------- Polynomial Model

% GenF18PolyModel; % bank angle can be changed here
% GenF18PolyBaselineModel;
% GenF18PolyRevisedModel; 

load 2_Mar14_PolyF18BaselineModel_Phi_35;

% Trim points
xtp = xtrim; utp = utrim; 

% ---- Generate Revised Cubic CLP Model
fcubicB = xcldotB;  
xclB = x; 



% ---- Generate Revised Cubic CLP Model
load 2_Mar14_PolyF18RevisedModel_Phi_35
fcubicR = xcldotR; %cleanpoly(xcldotR,[],1:3); 
xclR = x; 

% ---- initial condition
xx = [dxic(2:7);0]; 

% Simulate and check for convergence

event_params.nbig = 1e6;event_params.nsmall = 1e-6;
event_params.xbig = 1e3;event_params.xsmall = 1e-5;
event_params.funchandle = [];

% ----- Simulate Baseline Cubic System
[xtrajcubicB,xconv]= psim(fcubicB,xclB,xx,Tfinal,event_params);

% ----- Simulate Baseline Cubic System
[xtrajcubicR,xconv]= psim(fcubicR,xclR,xx,Tfinal,event_params);


% ---------------------------- Plot Results 

% ---------------------------- Baseline
figure5 = figure;
axes('Parent',figure5,'FontSize',18);
box('on');grid('on'); hold('all');


subplot1 = subplot(3,2,1,'FontSize',18);
hold on
plot(tB,yBresp(:,2)*r2d,'--b','LineWidth',3); 
plot(xtrajcubicB{1},(xtrajcubicB{2}(:,1)+xtp(2))*r2d,'-r','LineWidth',2);
ylabel('\beta (deg)','FontSize',18)
legend('Full','Cubic')
grid on 
box on 

subplot1 = subplot(3,2,2,'FontSize',18);
hold on
plot(tB,yBresp(:,4)*r2d,'--b','LineWidth',3);
plot(xtrajcubicB{1},(xtrajcubicB{2}(:,3)+ xtp(4))*r2d,'-r','LineWidth',2);
ylabel('p (deg/s)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,3,'FontSize',18);
hold on
plot(tB,yBresp(:,6)*r2d,'--b','LineWidth',3);
plot(xtrajcubicB{1},(xtrajcubicB{2}(:,5)+ xtp(6))*r2d,'-r','LineWidth',2);
ylabel('r (deg/s)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,4,'FontSize',18);
hold on
plot(tB,yBresp(:,7)*r2d,'--b','LineWidth',3);
plot(xtrajcubicB{1},(xtrajcubicB{2}(:,6)+ xtp(7))*r2d,'-r','LineWidth',2);
ylabel('\phi (deg)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,5,'FontSize',18);
hold on
plot(tB,yBresp(:,3)*r2d,'--b','LineWidth',3);
plot(xtrajcubicB{1},(xtrajcubicB{2}(:,2)+ xtp(3))*r2d,'-r','LineWidth',2);
ylabel('\alpha (deg)','FontSize',18)
xlabel('time, sec','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,6,'FontSize',18);
hold on
plot(tB,yBresp(:,5)*r2d,'--b','LineWidth',3);
plot(xtrajcubicB{1},(xtrajcubicB{2}(:,4)+ xtp(5))*r2d,'-r','LineWidth',2);
ylabel('q (deg/s)','FontSize',18)
xlabel('time, sec','FontSize',18)
grid on 
box on 

% ---------------------------- Revised
    
figure6 = figure;
axes('Parent',figure6,'FontSize',18);
box('on');grid('on'); hold('all');


subplot1 = subplot(3,2,1,'FontSize',18);
hold on
plot(tR,yRresp(:,2)*r2d,'--b','LineWidth',3); 
plot(xtrajcubicR{1},(xtrajcubicR{2}(:,1)+xtp(2))*r2d,'-r','LineWidth',2);
ylabel('\beta (deg)','FontSize',18)
legend('Full','Cubic')
grid on 
box on 

subplot1 = subplot(3,2,2,'FontSize',18);
hold on
plot(tR,yRresp(:,4)*r2d,'--b','LineWidth',3);
plot(xtrajcubicR{1},(xtrajcubicR{2}(:,3)+ xtp(4))*r2d,'-r','LineWidth',2);
ylabel('p (deg/s)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,3,'FontSize',18);
hold on
plot(tR,yRresp(:,6)*r2d,'--b','LineWidth',3);
plot(xtrajcubicR{1},(xtrajcubicR{2}(:,5)+ xtp(6))*r2d,'-r','LineWidth',2);
ylabel('r (deg/s)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,4,'FontSize',18);
hold on
plot(tR,yRresp(:,7)*r2d,'--b','LineWidth',3);
plot(xtrajcubicR{1},(xtrajcubicR{2}(:,6)+ xtp(7))*r2d,'-r','LineWidth',2);
ylabel('\phi (deg)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,5,'FontSize',18);
hold on
plot(tR,yRresp(:,3)*r2d,'--b','LineWidth',3);
plot(xtrajcubicR{1},(xtrajcubicR{2}(:,2)+ xtp(3))*r2d,'-r','LineWidth',2);
ylabel('\alpha (deg)','FontSize',18)
xlabel('time, sec','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,6,'FontSize',18);
hold on
plot(tR,yRresp(:,5)*r2d,'--b','LineWidth',3);
plot(xtrajcubicR{1},(xtrajcubicR{2}(:,4)+ xtp(5))*r2d,'-r','LineWidth',2);
ylabel('q (deg/s)','FontSize',18)
xlabel('time, sec','FontSize',18)
grid on 
box on 


