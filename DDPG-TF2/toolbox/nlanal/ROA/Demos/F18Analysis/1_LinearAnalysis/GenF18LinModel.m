%
% This file generates trim points for the F/A-18 aircraft. There could be
% three cases: 
%             (i)  Steady Straight & Level flight
%            (ii)  Steady Coordinated Turn 
%            (iii) UnCoordinated Turn 
%
% Moreover, this file also Generates the open loop 6-state polynomial model
%  around the trim point
%
% The initial condtions are described in x0. Refer to TrimF18.m file
% for more help.
% 
% Reference: (i) ETKIN,; Dynamics of Flight : Stability and Control 
%            (ii) A General Solution to the Aircraft Trim Problem
%
%
% Abhijit 11/17/2009

%==========================================================================
% Physical Data for F/A-18

 F18_AeroData = 1; 
 f18_data; 

%==========================================================================
% Cases for Trimming: 
% 
% Level = 0; 
% CoordTurn = 1;
% UnCoordTurn = 0; 


%==========================================================================
% Initial Condition 

if Level
    
    trimcase = 'level' ; 
    
    V0       =  450;         % Airspeed , ft/s
    beta0    =  0*d2r;       % Sideslip Angle, rad
    alpha0   =  3*d2r;       % Angle-of-attack, rad

    p0       =  0*d2r;       % Roll rate, rad/s
    q0       =  0*d2r;       % Pitch rate, rad/s
    r0       =  0*d2r;       % Yaw rate, rad/s

    phi0     =  0*d2r;       % Roll Angle, rad
    theta0   =  3*d2r;       % Pitch Angle, rad
    psi0     =  0*d2r;       % Yaw Angle, rad

    x0 = [V0; beta0; alpha0; p0; q0; r0; phi0; theta0; psi0]; 
    
elseif CoordTurn
    
    trimcase = 'coordturn' ;
    % Reference : A General Solution to the Aircraft Trim Problem
    
    V0       =  350;     % Airspeed , ft/s
    
    if ~exist('phi0','var')
        phi0     =  35*d2r;  % Roll Angle, rad
    end
    
    beta0    =  0*d2r;       % Sideslip Angle, rad
    alpha0   =  3*d2r;       % Angle-of-attack, rad
    
    p0       =  0*d2r;       % Roll rate, rad/s
    q0       =  0*d2r;       % Pitch rate, rad/s
    r0       =  0*d2r;       % Yaw rate, rad/s 
    theta0   =  3*d2r;       % Pitch Angle, rad
    psi0     =  0*d2r;       % Yaw Angle, rad
    
    x0 = [V0; beta0; alpha0; p0; q0; r0; phi0; theta0; psi0]; 
    
elseif UnCoordTurn

    trimcase = 'uncoordturn' ;
    % Reference : A General Solution to the Aircraft Trim Problem

    V0       =  350;     % Airspeed , ft/s
    
    if ~exist('phi0','var')
        phi0     =  35*d2r;  % Roll Angle, rad
    end
       
    beta0    =  10*d2r;      % Sideslip Angle, rad
    alpha0   =  3*d2r;       % Angle-of-attack, rad
    p0       =  0*d2r;       % Roll rate, rad/s
    q0       =  0*d2r;       % Pitch rate, rad/s
    r0       =  0*d2r;       % Yaw rate, rad/s
    theta0   =  3*d2r;       % Pitch Angle, rad
    psi0     =  0*d2r;       % Yaw Angle, rad

    x0 = [V0; beta0; alpha0; p0; q0; r0; phi0; theta0; psi0]; 

end


%==========================================================================
% Trim F/A-18 Aircraft 

displayon = 1;          % Display Trim values


% ------------- Set options for optimzation routine
% opt1 = optimset('MaxFunEvals',1e+04,'Display','off');
% opt = linoptions('OptimizationOptions',opt1,'DisplayReport','off');
opts.poly = 0; 

 if ~isfield(opts,'poly') || isempty(opts.poly)
       opts.poly    =  0; 
 end  
 
[xtrim, utrim, dx] = TrimF18(x0,opts, displayon,trimcase);


%==========================================================================
% Linearize F/A-18 Aircraft 

[A,B,C,D] = linmod('f18trim',xtrim,utrim);

% --------- Full 9-State Model
statenames = {'V (ft/s)','beta (rad)','alpha (rad)','p (rad/s)','q (rad/s)',...
            'r (rad/s)','phi (rad)','theta (rad)','psi (rad)'}';
inputnames = { 'aileron (rad)','rudder (rad)','stab (rad)', 'thrust (lbf)'}';       

linsys = ss(A,B,C,D,'statename',statenames,'inputname',inputnames,...
                    'outputname',statenames);

                
% --------- Reduced model (set thust to trim, ignore V, theta, psi) 
redstates = {'beta (rad)','alpha (rad)','p (rad/s)','q (rad/s)', ...
            'r (rad/s)','phi (rad)'}';

Ared = A(2:7,2:7); Bred = B(2:7,1:3); Cred = C(2:7,2:7); Dred = D(2:7,1:3); 
linsys_red = ss(Ared, Bred, Cred, Dred,'statename',redstates,...
                  'inputname',inputnames(1:3), 'outputname',redstates);
              
              
              
              
%==========================================================================
% Rationale of Decoupling V and theta.  Psi can be decoupled

if 0 

    % % ---- (i) See SVD Decomposition
    % [u0,s0,v0] = svd(A);  % unscaled
    % 
    % Scal = diag([15 5*d2r 10*d2r  20*d2r 10*d2r 5*d2r 15*d2r 10*d2r 5*d2r]); 
    % AScal = inv(Scal)*A*Scal; 
    % [u1,s1,v1] = svd(AScal);  % scaled

    % ---- (ii) bode analysis 

   
    wmin = 10^-2; 
    wmax = 10^2; 
    w = {wmin,wmax}; 
    w =logspace(wmin,wmax,25); 
  
    P = bodeoptions; 
    P.PhaseVisible = 'off';
    P.Grid = 'on'; 
    P.Title.String = '';  
    P.XLimmode = 'manual'; 
    P.XLim = [wmin wmax]; 
    P.YLimMode = 'auto'; 
    P.TickLabel.FontSize = 16; 
    P.XLabel.FontSize = 16; 
    P.YLabel.FontSize = 16; 
    P.OutputLabels.FontSize = 16;
    P.InputLabels.FontSize = 16; 
    P.Title.Color = [0 0 0]; 
    
    figure
    subplot(321)
    bodeplot(linsys(2,3),'-r',P)
    hold on 
    bodeplot(linsys_red(1,3),'--b',P)
    legend('9-State' , '6-State')
    
    subplot(322)
    bodeplot(linsys(3,3),'-r',P)
    hold on 
    bodeplot(linsys_red(2,3),'--b',P)

    
    subplot(323)
    P.Title.Color = [0 0 0];
    bodeplot(linsys(4,3),'-r',P)
    hold on 
    bodeplot(linsys_red(3,3),'--b',P)

    
    subplot(324)
    bodeplot(linsys(5,3),'-r',P)
    hold on 
    bodeplot(linsys_red(4,3),'--b',P)
       
    subplot(325)
    bodeplot(linsys(6,3),'-r',P)
    hold on 
    bodeplot(linsys_red(5,3),'--b',P)
  
    subplot(326)
    bodeplot(linsys(7,3),'-r',P)
    hold on 
    bodeplot(linsys_red(6,3),'--b',P)
    h = findobj(gcf,'type','line');
    set(h,'linewidth',2);
    
    h = get(gcf,'Children')
    % -------------------------------
    % Phase 
    
    P.PhaseVisible = 'on';
    P.MagVisible = 'off'; 
    P.Title.String = ''; 
    
    figure
    subplot(321)
    bodeplot(linsys(2,3),'-r',P)
    hold on 
    bodeplot(linsys_red(1,3),'--b',P)
    legend('9-State' , '6-State')
    
    subplot(322)
    bodeplot(linsys(3,3),'-r',P)
    hold on 
    bodeplot(linsys_red(2,3),'--b',P)

    
    subplot(323)
    bodeplot(linsys(4,3),'-r',P)
    hold on 
    bodeplot(linsys_red(3,3),'--b',P)
   
    subplot(324)
    bodeplot(linsys(5,3),'-r',P)
    hold on 
    bodeplot(linsys_red(4,3),'--b',P)
       
    subplot(325)
    bodeplot(linsys(6,3),'-r',P)
    hold on 
    bodeplot(linsys_red(5,3),'--b',P)
  
    subplot(326)
    bodeplot(linsys(7,3),'-r',P)
    hold on 
    bodeplot(linsys_red(6,3),'--b',P)
    h = findobj(gcf,'type','line');
    set(h,'linewidth',2);

   %---------------------------------------------------------------------
   % hand Drawn 
   
    wmin = 10^-2; 
    wmax = 10^2; 
    w0 = {wmin,wmax}; 
    % w =logspace(wmin,wmax,25); 
    [mag,pha,w] = bode(linsys(:,3),w0); 
    [mag_r,pha_r,w_r] = bode(linsys_red(:,3),w0); 
    
    figure
    subplot(321)
    semilogx(w,squeeze(20*log10(mag(2,:,:))),'-r')
    hold on 
    semilogx(w_r,squeeze(20*log10(mag_r(1,:,:))),'--b')
    grid on 
    set(gca,'FontSize',18)
    ylabel('Magnitude (dB)','FontSize',18)
    title('\delta_{stab}(rad) to \beta (rad)', 'FontSize',16)
    legend('9-State' , '6-State')
    
    subplot(322)
    semilogx(w,squeeze(20*log10(mag(3,:,:))),'-r')
    hold on 
    semilogx(w_r,squeeze(20*log10(mag_r(2,:,:))),'--b')
    set(gca,'FontSize',18)
    ylabel('Magnitude (dB)','FontSize',18)
    title('\delta_{stab}(rad) to \alpha (rad)', 'FontSize',16)
    grid on 

    subplot(323)
    semilogx(w,squeeze(20*log10(mag(4,:,:))),'-r')
    hold on 
    semilogx(w_r,squeeze(20*log10(mag_r(3,:,:))),'--b')
    set(gca,'FontSize',18)
    ylabel('Magnitude (dB)','FontSize',18)
    title('\delta_{stab}(rad) to p (rad/s)', 'FontSize',16)
    grid on 
 
    subplot(324)
    semilogx(w,squeeze(20*log10(mag(5,:,:))),'-r')
    hold on 
    semilogx(w_r,squeeze(20*log10(mag_r(4,:,:))),'--b')
    set(gca,'FontSize',18)
    ylabel('Magnitude (dB)','FontSize',18)
    title('\delta_{stab}(rad) to q (rad/s)', 'FontSize',16)
    grid on 
 
    subplot(325)
    semilogx(w,squeeze(20*log10(mag(6,:,:))),'-r')
    hold on 
    semilogx(w_r,squeeze(20*log10(mag_r(5,:,:))),'--b')
    set(gca,'FontSize',18)
    ylabel('Magnitude (dB)','FontSize',18)
    xlabel('Frequency (rad/s)','FontSize',18)
    title('\delta_{stab}(rad) to r (rad/s)', 'FontSize',16)
    grid on 

    subplot(326)
    semilogx(w,squeeze(20*log10(mag(7,:,:))),'-r')
    hold on 
    semilogx(w_r,squeeze(20*log10(mag_r(6,:,:))),'--b')
    set(gca,'FontSize',18)
    xlabel('Frequency (rad/s)','FontSize',18)
    ylabel('Magnitude (dB)','FontSize',18)
    title('\delta_{stab}(rad) to \phi (rad)', 'FontSize',16)
    grid on 
    
    h = findobj(gcf,'type','line');
    set(h,'linewidth',2);
    
   
    figure
    subplot(321)
    semilogx(w,squeeze(pha(2,:,:)),'-r')
    hold on 
    semilogx(w_r,squeeze(pha_r(1,:,:)),'--b')
    grid on 
    set(gca,'FontSize',18)
    ylabel('Phase (deg)','FontSize',18)
    title('\delta_{stab}(rad) to \beta (rad)', 'FontSize',16)
    legend('9-State' , '6-State')
    
    subplot(322)
    semilogx(w,squeeze(pha(3,:,:)),'-r')
    hold on 
    semilogx(w_r,squeeze(pha_r(2,:,:)),'--b')
    set(gca,'FontSize',18)
    ylabel('Phase (deg)','FontSize',18)
    title('\delta_{stab}(rad) to \alpha (rad)', 'FontSize',16)
    grid on 

    subplot(323)
    semilogx(w,squeeze(pha(4,:,:)),'-r')
    hold on 
    semilogx(w_r,squeeze(pha_r(3,:,:)),'--b')
    set(gca,'FontSize',18)
    ylabel('Phase (deg)','FontSize',18)
    title('\delta_{stab}(rad) to p (rad/s)', 'FontSize',16)
    grid on 
 
    subplot(324)
    semilogx(w,squeeze(pha(5,:,:)),'-r')
    hold on 
    semilogx(w_r,squeeze(pha_r(4,:,:)),'--b')
    set(gca,'FontSize',18)
    ylabel('Phase (deg)','FontSize',18)
    title('\delta_{stab}(rad) to q (rad/s)', 'FontSize',16)
    grid on 
 
    subplot(325)
    semilogx(w,squeeze(pha(6,:,:)),'-r')
    hold on 
    semilogx(w_r,squeeze(pha_r(5,:,:)),'--b')
    set(gca,'FontSize',18)
    ylabel('Phase (deg)','FontSize',18)
    xlabel('Frequency (rad/s)','FontSize',18)
    title('\delta_{stab}(rad) to r (rad/s)', 'FontSize',16)
    grid on 

    subplot(326)
    semilogx(w,squeeze(pha(7,:,:)),'-r')
    hold on 
    semilogx(w_r,squeeze(pha_r(6,:,:)),'--b')
    set(gca,'FontSize',18)
    xlabel('Frequency (rad/s)','FontSize',18)
    ylabel('Phase (deg)','FontSize',18)
    title('\delta_{stab}(rad) to \phi (rad)', 'FontSize',16)
    grid on 
    
    h = findobj(gcf,'type','line');
    set(h,'linewidth',2);
    

end