%
% This file compares the Closed Loop Models for initial condition perturbation  
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


% ---------------------------- Plot Results 

% ---------------------------- Baseline
figure5 = figure;
axes('Parent',figure5,'FontSize',18);
box('on');grid('on'); hold('all');


subplot1 = subplot(3,2,1,'FontSize',18);
hold on
plot(tB,yBresp(:,2)*r2d,'--r','LineWidth',3); 
plot(tR,yRresp(:,2)*r2d,'-b','LineWidth',3); 
ylabel('\beta (deg)','FontSize',18)
legend('Baseline','Revised')
grid on 
box on 

subplot1 = subplot(3,2,2,'FontSize',18);
hold on
plot(tB,yBresp(:,4)*r2d,'--r','LineWidth',3);
plot(tR,yRresp(:,4)*r2d,'-b','LineWidth',3);
ylabel('p (deg/s)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,3,'FontSize',18);
hold on
plot(tB,yBresp(:,6)*r2d,'--r','LineWidth',3);
plot(tR,yRresp(:,6)*r2d,'-b','LineWidth',3);
ylabel('r (deg/s)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,4,'FontSize',18);
hold on
plot(tB,yBresp(:,7)*r2d,'--r','LineWidth',3);
plot(tR,yRresp(:,7)*r2d,'-b','LineWidth',3);
ylabel('\phi (deg)','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,5,'FontSize',18);
hold on
plot(tB,yBresp(:,3)*r2d,'--r','LineWidth',3);
plot(tR,yRresp(:,3)*r2d,'-b','LineWidth',3);
ylabel('\alpha (deg)','FontSize',18)
xlabel('time, sec','FontSize',18)
grid on 
box on 

subplot1 = subplot(3,2,6,'FontSize',18);
hold on
plot(tB,yBresp(:,5)*r2d,'--r','LineWidth',3);
plot(tB,yBresp(:,5)*r2d,'-b','LineWidth',3);
ylabel('q (deg/s)','FontSize',18)
xlabel('time, sec','FontSize',18)
grid on 
box on 

