%
% This file Simulate F/A-18 Models. The following simulations are performed: 
% (1) Falling Leaf Generation : F18 Full 9-State Model  
%
% Abhijit December 09

 

% F18 Data
f18_data

if 0 

% ===================================================================
% (1) Falling Leaf Generation : F18 Full 9-State Model


% State Initial Value
% V       =  350;          % Airspeed , ft/s
% beta    =  20*d2r;       % Sideslip Angle, rad
% alpha   =  40*d2r;       % Angle-of-attack, rad
% 
% p       =  10*d2r;       % Roll rate, rad/s
% q       =  0*d2r;       % Pitch rate, rad/s
% r       =  5*d2r;       % Yaw rate, rad/s
% 
% phi     =  0*d2r;       % Roll Angle, rad
% theta   =  0*d2r;       % Pitch Angle, rad
% psi     =  0*d2r;       % Yaw Angle, rad

V       =  350;          % Airspeed , ft/s
beta    =  20*d2r;       % Sideslip Angle, rad
alpha   =  40*d2r;       % Angle-of-attack, rad

p       =  10*d2r;       % Roll rate, rad/s
q       =  0*d2r;       % Pitch rate, rad/s
r       =  5*d2r;       % Yaw rate, rad/s

phi     =  0*d2r;       % Roll Angle, rad
theta   =  0*d2r;       % Pitch Angle, rad
psi     =  0*d2r;       % Yaw Angle, rad

% Stack Initial Condition for State
xfl = [V;beta;alpha;p;q;r;phi;theta;psi];  

tfl = [0 25];                % time scale 
ufl = [0;zeros(3,1)];        % Input Initialization

%---------------------------------------------
% Simulate
opts.poly = 0; 
[t,y] = ode45(@f18full,tfl,xfl,[],ufl,opts); 


%-------------------------------------------------
% Plot Results

figure1 = figure;
axes('Parent',figure1,'FontSize',14);
plot(t,y(:,4)*r2d,'-r','LineWidth',2)
hold on 
plot(t,y(:,6)*r2d,'--g','LineWidth',2)
grid on 
legend('Roll Rate (deg/s)','Yaw Rate (deg/s)')
xlabel('time, sec')
ylabel('Rate (deg/s)')


figure2 = figure;
axes('Parent',figure2,'FontSize',14);
plot(y(:,2)*r2d,y(:,3)*r2d)
xlabel('sideslip (deg)')
ylabel('aoa (deg)')
grid on 


% ===================================================================
% (2) Compare Simulation Comparison with Baseline and Revised

Bank25 = load('F18Model_Coord_Phi_35','Pss','xtrim','utrim','CssB','CssR'); 

% ----- Extract information needed for simulink model
xtrim = Bank25.xtrim; 
utrim = Bank25.utrim; 
C = Bank25.Pss.C; 
D = Bank25.Pss.D; 
opts.poly = 0;

% ------- Sim paramter setup
aytrim = C(1,:)*xtrim(2:7) + D(1,:)*utrim(1:3);
bdottrim = C(end,:)*xtrim(2:7) + D(end,:)*utrim(1:3);


InpPert = 1;   % Perturb the input channel 
dx0 = xtrim;   % start from the trim state

% ---- Baseline Simulation 

K = Bank25.CssB; 
[tB, yB] = sim('f18nlclpsim',[0 10]); 

K = Bank25.CssR; 
[tR, yR] = sim('f18nlclpsim',[0 10]); 

% ------ Plot Results 

figure

subplot(321)
ph1 = plot(tB,yB(:,2)*r2d,'--r',tR,yR(:,2)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('beta (deg)','FontSize',14)
set(gca,'FontSize',14); 
grid on 

subplot(323)
ph1 = plot(tB,yB(:,3)*r2d,'--r',tR,yR(:,3)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('alpha (deg)','FontSize',14)
set(gca,'FontSize',14); grid on 

subplot(325)
ph1 = plot(tB,yB(:,7)*r2d,'--r',tR,yR(:,7)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('phi (deg)','FontSize',14)
set(gca,'FontSize',14); grid on 

subplot(322)
ph1 = plot(tB,yB(:,4)*r2d,'--r',tR,yR(:,4)*r2d,'-b'); 
set(ph1,'LineWidth',2)
ylabel('p (deg/s)','FontSize',14)
set(gca,'FontSize',14); grid on 
lh = legend('Baseline' , ' Revised'); 
set(lh, 'FontSize',14)


% subplot(322)
% ph1 = plot(tB,yB(:,4)*r2d,'--r'); 
% hold on 
% [AX,H1,H2]=plotyy(tR,yR(:,4)*r2d,tR,inp*r2d);
% set(ph1,'LineWidth',2)
% set(get(AX(2),'Ylabel'),'String','ail_{cmd} (deg)') 
% set(get(AX(1),'Ylabel'),'String','p (deg/s)') 
% set(H1,'LineStyle','-','LineWidth',2)
% set(H2,'LineStyle',':','LineWidth',2)
% grid on 
% lh = legend('Baseline' , ' Revised', 'Commanded Input');
% set(lh, 'FontSize',14)

subplot(324)
ph1 = plot(tB,yB(:,5)*r2d,'--r',tR,yR(:,5)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('q (deg/s)','FontSize',14)
set(gca,'FontSize',14); grid on 

subplot(326)
ph1 = plot(tB,yB(:,6)*r2d,'--r',tR,yR(:,6)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('r (deg/s)','FontSize',14)
set(gca,'FontSize',14); grid on 

figure
% subplot(311)
% ph1 = plot(tB,yB(:,10)*r2d,'--r',tR,yR(:,10)*r2d,'-b');
% hold on 
% plot(tR,inp*r2d,'.-r')
% set(ph1,'LineWidth',2)
% ylabel('ail_{act} (deg)','FontSize',14)
% set(gca,'FontSize',14); 
% grid on 
% lh = legend('Baseline' , ' Revised', 'Commanded Input');
% set(lh, 'FontSize',14)

subplot(311)

ph1 = plot(tB,yB(:,10)*r2d,'--r'); 
hold on 
[AX,H1,H2]=plotyy(tR,yR(:,10)*r2d,tR,inp*r2d);
set(ph1,'LineWidth',2)
set(get(AX(2),'Ylabel'),'String','ail_{cmd} (deg)') 
set(get(AX(1),'Ylabel'),'String','ail_{act} (deg)') 
set(H1,'LineStyle','-','LineWidth',2)
set(H2,'LineStyle',':','LineWidth',4)
grid on 
axes(AX(1)); set(gca,'fontsize',18);
axes(AX(2)); set(gca,'fontsize',18);
lh = legend('Commanded Input','Baseline' , ' Revised');
set(lh, 'FontSize',14)


subplot(312)
ph1 = plot(tB,yB(:,11)*r2d,'--r',tR,yR(:,11)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('rud_{act} (deg)','FontSize',14)
set(gca,'FontSize',14); grid on 

subplot(313)
ph1 = plot(tB,yB(:,12)*r2d,'--r',tR,yR(:,12)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('stab_{act} (deg)','FontSize',14)
set(gca,'FontSize',14); grid on 




else

%==========================================================================
% Simulate Worst Case perturbation in nonlinear sims

Bank25 = load('F18Model_Coord_Phi_35','Pss','xtrim','utrim','CssB','CssR'); 

% ----- Extract information needed for simulink model
xtrim = Bank25.xtrim; 
utrim = Bank25.utrim; 
C = Bank25.Pss.C; 
D = Bank25.Pss.D; 
opts.poly = 0;

% ------- Sim paramter setup
aytrim = C(1,:)*xtrim(2:7) + D(1,:)*utrim(1:3);
bdottrim = C(end,:)*xtrim(2:7) + D(end,:)*utrim(1:3);


InpPert = 0;   % Perturb the input channel 
dx0 = xtrim;   % start from the trim state

% ---- Run 'RunF18_WCSim.m' 

% ------ Baseline WC Perturbation
opts.d1= 0.0299 -  Bank25.Pss.A(1,1);  
opts.d2= -6.4608 - Bank25.Pss.A(3,1);  
opts.d3= 0.3836  - Bank25.Pss.A(5,1);  
opts.d4= -0.1440 - Bank25.Pss.A(5,5);
opts.d5= -0.5555 - Bank25.Pss.A(3,3); 
opts.d6 = 1.1000 - Bank25.Pss.A(2,4); 
opts.d7= -0.7801 - Bank25.Pss.A(4,2); 
opts.d8= -0.2142 - Bank25.Pss.A(4,4);

% opts.d1= 0;  
% opts.d2= 0;  
% opts.d3= 0;  
% opts.d4= 0; 
% opts.d5= 0; 
% opts.d6 = 0;   
% opts.d7= 0;  
% opts.d8= 0;


load WCGainIMUDiag; 
d11 = maxgainuncBBeta.d11;  d22 = maxgainuncBBeta.d22; d33 = maxgainuncBBeta.d33; 
d11 = ss(0) ; d22 = ss(0); d33 = ss(0); 

% ---- Baseline Simulation 

K = Bank25.CssB; 
[tB, yB] = sim('f18nlclpWCsim',[0 25]); 


% ------ Revised WC Perturbation
opts.d1= 0.0299 - Bank25.Pss.A(1,1);
opts.d2= -6.4608 - Bank25.Pss.A(3,1); 
opts.d3= 0.3836 - Bank25.Pss.A(5,1); 
opts.d4= -0.1440 - Bank25.Pss.A(5,5);
opts.d5= -0.6790 - Bank25.Pss.A(3,3); 
opts.d6= 0.9000 - Bank25.Pss.A(2,4);
opts.d7= -0.9534 - Bank25.Pss.A(4,2); 
opts.d8= -0.1752 - Bank25.Pss.A(4,4);

% opts.d1= 0;  
% opts.d2= 0;  
% opts.d3= 0;  
% opts.d4= 0; 
% opts.d5= 0; 
% opts.d6 = 0;   
% opts.d7= 0;  
% opts.d8= 0;

d11 = maxgainuncRBeta.d11;  d22 = maxgainuncRBeta.d22;  d33 = maxgainuncRBeta.d33; 
d11 = ss(0) ; d22 = ss(0); d33 = ss(0); 

K = Bank25.CssR; 
[tR, yR] = sim('f18nlclpWCsim',[0 25]); 

% ------ Plot Results 

figure

subplot(321)
ph1 = plot(tB,yB(:,2)*r2d,'--r',tR,yR(:,2)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('beta (deg)','FontSize',14)
set(gca,'FontSize',18); 
grid on 

subplot(323)
ph1 = plot(tB,yB(:,3)*r2d,'--r',tR,yR(:,3)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('alpha (deg)','FontSize',14)
set(gca,'FontSize',18); grid on 

subplot(325)
ph1 = plot(tB,yB(:,7)*r2d,'--r',tR,yR(:,7)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('phi (deg)','FontSize',14)
set(gca,'FontSize',18); grid on 

subplot(322)
ph1 = plot(tB,yB(:,4)*r2d,'--r',tR,yR(:,4)*r2d,'-b'); 
set(ph1,'LineWidth',2)
ylabel('p (deg/s)','FontSize',14)
set(gca,'FontSize',14); grid on 
lh = legend('Baseline' , ' Revised'); 
set(lh, 'FontSize',18)

subplot(324)
ph1 = plot(tB,yB(:,5)*r2d,'--r',tR,yR(:,5)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('q (deg/s)','FontSize',14)
set(gca,'FontSize',18); grid on 

subplot(326)
ph1 = plot(tB,yB(:,6)*r2d,'--r',tR,yR(:,6)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('r (deg/s)','FontSize',14)
set(gca,'FontSize',18); grid on 


figure
subplot(311)
ph1 = plot(tB,yB(:,10)*r2d,'--r',tR,yR(:,10)*r2d,'-b'); 
set(ph1,'LineWidth',2)
ylabel('ail_{act} (deg)','FontSize',14)
set(gca,'FontSize',18); grid on 



subplot(312)
ph1 = plot(tB,yB(:,11)*r2d,'--r',tR,yR(:,11)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('rud_{act} (deg)','FontSize',14)
set(gca,'FontSize',18); grid on 

subplot(313)
ph1 = plot(tB,yB(:,12)*r2d,'--r',tR,yR(:,12)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('stab_{act} (deg)','FontSize',14)
set(gca,'FontSize',18); grid on 


figure
ph1 = plot(tB,yB(:,4)*r2d,'--r',tB,yB(:,6)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('p-r (deg)','FontSize',14)
set(gca,'FontSize',14); 
grid on 

figure
ph1 = plot(yB(:,2)*r2d,yB(:,3)*r2d,'-b');
set(ph1,'LineWidth',2)
ylabel('beta-aoa (deg)','FontSize',14)
set(gca,'FontSize',14); 
grid on 



end