% ========================================================================
%  Monte Carlo Simulation Plot

d2r = pi/180; r2d = 1 / d2r; 
Dmax = diag([10*d2r 25*d2r 35*d2r 30*d2r 15*d2r 25*d2r 20*d2r]);  
N = inv(Dmax^2);   %  Scale by inverse of max state values

% ------------------------------------------------------------
% Baseline 
if 0 
    
load 2_Mar14_PolyF18BaselineModel_Phi_35
f = cleanpoly(xcldotB,10^-06);  

BClawMC = load('1_Mar14_ROAMC_F18Baseline_35'); 
[b , i] = min(BClawMC.betaun);
x01 = BClawMC.xun(:,i); 
b1 = x01'*N*x01 ;
scl = sqrt( (x01'*N*x01)/0.985/b1 );
x01 = x01/scl;
b1= x01'*N*x01

Tfinal = 30; 
event_params.nbig = 1e6; event_params.nsmall = 1e-6;
event_params.xbig = 1000; event_params.xsmall = 1e-5;
event_params.funchandle = [];

[xtrajBu,xconv]= psim(f,x,x01,Tfinal,event_params);
dailBu = 0.08*xtrajBu{2}(:,3);
drudBu = -xtrajBu{2}(:,7)  - 1.1*xtrajBu{2}(:,5) - ...
    0.5.*(-0.15101*xtrajBu{2}(:,2).^2.*xtrajBu{2}(:,1) + ...
    0.20808*xtrajBu{2}(:,2).*xtrajBu{2}(:,1) - 0.5758*xtrajBu{2}(:,1)); 
dstabBu = -0.8*xtrajBu{2}(:,2) - 8*xtrajBu{2}(:,4); 



    
figure5 = figure;
axes('Parent',figure5,'FontSize',16);
box('on');grid('on'); hold('all');

subplot1 = subplot(3,1,1,'FontSize',16);
hold on
plot(xtrajBu{1}, xtrajBu{2}(:,1)*r2d,'--b','LineWidth',2)
plot(xtrajBu{1}, xtrajBu{2}(:,2)*r2d,'-b','LineWidth',2)
plot(xtrajBu{1}, xtrajBu{2}(:,6)*r2d,'-.b','LineWidth',2)
legend('\beta', '\alpha', '\phi','Fontsize',16)
ylabel('State (deg)','Fontsize',16)
grid on 

subplot1 = subplot(3,1,2,'FontSize',16);
hold on
plot(xtrajBu{1}, xtrajBu{2}(:,3)*r2d,'--b','LineWidth',2)
plot(xtrajBu{1}, xtrajBu{2}(:,4)*r2d,'-b','LineWidth',2)
plot(xtrajBu{1}, xtrajBu{2}(:,5)*r2d,'-.b','LineWidth',2)
legend('p', 'q', 'r','Fontsize',16)
ylabel(' Rates (deg/s)','Fontsize',16)
grid on 

subplot1 = subplot(3,1,3,'FontSize',16);
hold on
plot(xtrajBu{1}, -dailBu*r2d,'--b','LineWidth',2)
plot(xtrajBu{1}, -drudBu*r2d,'-b','LineWidth',2)
plot(xtrajBu{1}, -dstabBu*r2d,'-.b','LineWidth',2)
legend('\delta_{ail}', '\delta_{rud}', '\delta_{stab}','Fontsize',16)
ylabel('Control Signal (deg)','Fontsize',16)
xlabel('time (sec)','Fontsize',16)
grid on


% -------- Stable baseline Scaling 
Tfinal = 30; 
scl = sqrt( (x01'*N*x01)/0.995/b1 );
x0 = x01/scl;
b2= x0'*N*x0

[xtrajBs,xconv]= psim(f,x,x0,Tfinal,event_params);
dailBs = 0.08*xtrajBs{2}(:,3);
drudBs = -xtrajBs{2}(:,7)  - 1.1*xtrajBs{2}(:,5) - ...
    0.5.*(-0.15101*xtrajBs{2}(:,2).^2.*xtrajBs{2}(:,1) + ...
    0.20808*xtrajBs{2}(:,2).*xtrajBs{2}(:,1) - 0.5758*xtrajBs{2}(:,1)); 
dstabBs = -0.8*xtrajBs{2}(:,2) - 8*xtrajBs{2}(:,4); 

figure5 = figure;
axes('Parent',figure5,'FontSize',16);
box('on');grid('on'); hold('all');

subplot1 = subplot(3,1,1,'FontSize',16);
hold on
plot(xtrajBs{1}, xtrajBs{2}(:,1)*r2d,'--b','LineWidth',2)
plot(xtrajBs{1}, xtrajBs{2}(:,2)*r2d,'-b','LineWidth',2)
plot(xtrajBs{1}, xtrajBs{2}(:,6)*r2d,'-.b','LineWidth',2)
legend('\beta', '\alpha', '\phi','Fontsize',16)
ylabel('State (deg)','Fontsize',16)
grid on 

subplot1 = subplot(3,1,2,'FontSize',16);
hold on
plot(xtrajBs{1}, xtrajBs{2}(:,3)*r2d,'--b','LineWidth',2)
plot(xtrajBs{1}, xtrajBs{2}(:,4)*r2d,'-b','LineWidth',2)
plot(xtrajBs{1}, xtrajBs{2}(:,5)*r2d,'-.b','LineWidth',2)
legend('p', 'q', 'r','Fontsize',16)
ylabel(' Rates (deg/s)','Fontsize',16)
grid on 

subplot1 = subplot(3,1,3,'FontSize',16);
hold on
plot(xtrajBs{1}, -dailBs*r2d,'--b','LineWidth',2)
plot(xtrajBs{1}, -drudBs*r2d,'-b','LineWidth',2)
plot(xtrajBs{1}, -dstabBs*r2d,'-.b','LineWidth',2)
legend('\delta_{ail}', '\delta_{rud}', '\delta_{stab}','Fontsize',16)
ylabel('Control Signal (deg)','Fontsize',16)
xlabel('time (sec)','Fontsize',16)
grid on


else
% ------------------------------------------------------------
% Revised

load 2_Mar14_PolyF18RevisedModel_Phi_35
f = cleanpoly(xcldotR,10^-06);
 
RClawMC = load('1_Mar14_ROAMC_F18Revised_35'); 
[b , i] = min(RClawMC.betaun);
x01  = RClawMC.xun(:,i); 
b1 = x01'*N*x01 ;
% scl = sqrt( (x01'*N*x01)/0.99/b1 );
% x01 = x01/scl;
% b1= x01'*N*x01

Tfinal = 10; 
event_params.nbig = 1e6; event_params.nsmall = 1e-6;
event_params.xbig = 1000; event_params.xsmall = 1e-5;
event_params.funchandle = [];

[xtrajRu,xconv]= psim(f,x,x01,Tfinal,event_params);
dailRu = 0.08*xtrajRu{2}(:,3);
drudRu = -xtrajRu{2}(:,7)  - 1.1*xtrajRu{2}(:,5) - ...
    0.5.*(-0.15101*xtrajRu{2}(:,2).^2.*xtrajRu{2}(:,1) + ...
    0.20808*xtrajRu{2}(:,2).*xtrajRu{2}(:,1) - 0.5758*xtrajRu{2}(:,1)); 
dstabRu = -0.8*xtrajRu{2}(:,2) - 8*xtrajRu{2}(:,4); 



figure5 = figure;
axes('Parent',figure5,'FontSize',16);
box('on');grid('on'); hold('all');

subplot1 = subplot(3,1,1,'FontSize',16);
hold on
plot(xtrajRu{1}, xtrajRu{2}(:,1)*r2d,'--b','LineWidth',2)
plot(xtrajRu{1}, xtrajRu{2}(:,2)*r2d,'-b','LineWidth',2)
plot(xtrajRu{1}, xtrajRu{2}(:,6)*r2d,'-.b','LineWidth',2)
legend('\beta', '\alpha', '\phi','Fontsize',16)
ylabel('State (deg)','Fontsize',16)
grid on 

subplot1 = subplot(3,1,2,'FontSize',16);
hold on
plot(xtrajRu{1}, xtrajRu{2}(:,3)*r2d,'--b','LineWidth',2)
plot(xtrajRu{1}, xtrajRu{2}(:,4)*r2d,'-b','LineWidth',2)
plot(xtrajRu{1}, xtrajRu{2}(:,5)*r2d,'-.b','LineWidth',2)
legend('p', 'q', 'r','Fontsize',16)
ylabel(' Rates (deg/s)','Fontsize',16)
grid on 

subplot1 = subplot(3,1,3,'FontSize',16);
hold on
plot(xtrajRu{1}, -dailRu*r2d,'--b','LineWidth',2)
plot(xtrajRu{1}, -drudRu*r2d,'-b','LineWidth',2)
plot(xtrajRu{1}, -dstabRu*r2d,'-.b','LineWidth',2)
legend('\delta_{ail}', '\delta_{rud}', '\delta_{stab}','Fontsize',16)
ylabel('Control Signal (deg)','Fontsize',16)
xlabel('time (sec)','Fontsize',16)
grid on



% -------- Stable revised Scaling 
Tfinal = 5; 
scl = sqrt( (x01'*N*x01)/0.995/b1 );
x0 = x01/scl;
x0'*N*x0 

[xtrajRs,xconv]= psim(f,x,x0,Tfinal,event_params);
dailRs = 0.08*xtrajRs{2}(:,3);
drudRs = -xtrajRs{2}(:,7)  - 1.1*xtrajRs{2}(:,5) - ...
    0.5.*(-0.15101*xtrajRs{2}(:,2).^2.*xtrajRs{2}(:,1) + ...
    0.20808*xtrajRs{2}(:,2).*xtrajRs{2}(:,1) - 0.5758*xtrajRs{2}(:,1)); 
dstabRs = -0.8*xtrajRs{2}(:,2) - 8*xtrajRs{2}(:,4); 

figure5 = figure;
axes('Parent',figure5,'FontSize',16);
box('on');grid('on'); hold('all');

subplot1 = subplot(3,1,1,'FontSize',16);
hold on
plot(xtrajRs{1}, xtrajRs{2}(:,1)*r2d,'--b','LineWidth',2)
plot(xtrajRs{1}, xtrajRs{2}(:,2)*r2d,'-b','LineWidth',2)
plot(xtrajRs{1}, xtrajRs{2}(:,6)*r2d,'-.b','LineWidth',2)
legend('\beta', '\alpha', '\phi','Fontsize',16)
ylabel('State (deg)','Fontsize',16)
grid on 

subplot1 = subplot(3,1,2,'FontSize',16);
hold on
plot(xtrajRs{1}, xtrajRs{2}(:,3)*r2d,'--b','LineWidth',2)
plot(xtrajRs{1}, xtrajRs{2}(:,4)*r2d,'-b','LineWidth',2)
plot(xtrajRs{1}, xtrajRs{2}(:,5)*r2d,'-.b','LineWidth',2)
legend('p', 'q', 'r','Fontsize',16)
ylabel(' Rates (deg/s)','Fontsize',16)
grid on 

subplot1 = subplot(3,1,3,'FontSize',16);
hold on
plot(xtrajRs{1}, -dailRs*r2d,'--b','LineWidth',2)
plot(xtrajRs{1}, -drudRs*r2d,'-b','LineWidth',2)
plot(xtrajRs{1}, -dstabRs*r2d,'-.b','LineWidth',2)
legend('\delta_{ail}', '\delta_{rud}', '\delta_{stab}','Fontsize',16)
ylabel('Control Signal (deg)','Fontsize',16)
xlabel('time (sec)','Fontsize',16)
grid on

end