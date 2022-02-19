%
% This file generates the Open-Loop / Closed-loop linear model for both 
% baseline and revised model


% ------------------------------------------------------------
% Form the appropiate output and input

pvar beta alpha p q r phi d_AIL d_RUD d_STAB 

% ------------------ State Derivative / Input / Output / States
x    = [beta;  alpha;   p;   q;    r;    phi];           % State
u    = [d_AIL;  d_RUD;  d_STAB];                         % Input


F18_AeroData = 1; f18_data;

Cyb     = F18Aero.Cyb_0 + F18Aero.Cyb_1*alpha + F18Aero.Cyb_2*alpha^2 ; 
Cyda    = F18Aero.Cyda_0 + F18Aero.Cyda_1*alpha + F18Aero.Cyda_2*alpha^2 ...
         + F18Aero.Cyda_3*alpha^3; 
Cydr    = F18Aero.Cydr_0 + F18Aero.Cydr_1*alpha + F18Aero.Cydr_2*alpha^2 ...
         + F18Aero.Cydr_3*alpha^3; 

% -------- Total Side Force
C_Y     = Cyb*beta + Cydr*d_RUD +  Cyda*d_AIL;
  
% ---------- Lateral Acceleration
Vtrim = xtrim(1); 
qbar  =  (1/2)*rho*Vtrim^2; 
Ay = qbar*S*(C_Y)/m/g; 

% ---------- Linear Estimation of betadot term
betalind = Ared(1,:)*x + Bred(1,:)*u;

% ------ Output Suitable for both baseline and revised model
y = [Ay; p; r; alpha; beta; q; betalind];

% ------- Form open Loop Linear SS
A0 = Ared; B0 = Bred; 
[C0 , D0 ] = plinearize( y,x,u);

statenames = {'beta (rad)','alpha (rad)','p (rad/s)','q (rad/s)',...
                'r (rad/s)','phi (rad)'}';
inputnames = {'aileron (rad)','rudder (rad)','stab (rad)'}';       
outputnames = {'a_y (g)','p (rad/s)','r (rad/s)','alpha (rad)', ...
                'beta (rad)','q (rad/s)','betadot'}';        
Pss = ss(A0, B0, C0, D0, 'statename',statenames,...
        'inputname',inputnames,'outputname',outputnames);
    
    
% =======================================================================
%  Form Baseline Controller
%
%  xc = [xcB];
%  u = [d_AIL; d_RUD; d_STAB];
%  y = [ay; p; r; alpha;0; q; 0];
%
%  xc_dot = Ac*xc+Bc*y
%  u      = Cc*xc+Dc*y


%------------- Baseline Controller Realization
AcB = [ -1 ];
BcB = -[0 0 -4.9 0 0 0 0]; 
CcB = [0;-1;0];
DcB = -[0 -0.08 0 0 0 0 0; 0.5 0 1.1 0 0 0 0;  0 0 0 0.8 0 8 0];
CssB = ss(AcB,BcB,CcB,DcB);


% =======================================================================
%  Form Revised Controller
%
%  xc = [xcR];
%  u =  [d_AIL; d_RUD; d_STAB];
%  y =  [ay; p; r; alpha; beta; q; betalind];
%
%  xc_dot = Ac*xc+Bc*y
%  u      = Cc*xc+Dc*u

kbetadot = 2;              % beta dot feedback
kbeta    = 0.5;            % beta feedback

% Revised Controller Realization
AcR = [-1 ];
BcR = -[0 0 -4.9 0 0 0 0]; 
CcR =  [0 ;-1; 0];
DcR = -[0 -0.08 0 0 -kbeta 0 -kbetadot; 0.5 0 1.1 0 0 0 0;  0 0 0 0.8 0 8 0];
CssR = ss(AcR,BcR,CcR,DcR);

%==========================================================================
% Linearized Closed Loop  Plant

% --------- Baseline : Negative Feedback of Pss.Baseline and CssB
CLB = feedback(Pss,CssB);

% --------- Baseline : Negative Feedback of Pss.Baseline and CssB
CLR = feedback(Pss,CssR);

fprintf('\n ======================================= \n')
fprintf('Baseline with phi = %d  \t beta = %4.6f\n',xtrim(7)*r2d,xtrim(2)*r2d)
fprintf('\n \n')
damp(CLB)

fprintf('\n ======================================= \n')
fprintf('Revised with phi = %d  \t beta = %4.6f\n',xtrim(7)*r2d,xtrim(2)*r2d)
fprintf('\n \n')
damp(CLR)