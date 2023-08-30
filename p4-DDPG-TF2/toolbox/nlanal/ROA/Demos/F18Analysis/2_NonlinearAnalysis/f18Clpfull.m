function [xcldot2 , xcldot, F, G] = f18Clpfull(t,x,opts)
%
% [xcldot2 , xcldot, F, G] = f18Clpfull(t,x,opts)
%
% This function creates a closed-loop realization of the F18 plant around a
% trim point. this realization does not involve any approximation, it is
% based on full nonlinear model. The F/A-18 plant has been reduced to six
% state representation of open loop as opposed to 9-state. The plant is 
% described as 
%               xdot = F(x)   + G(x)*u 
%               y    = H(x)   + J(x)*u
%
% Controllers is described as: 
%               xcdot = Ac*xc + Bc*y
%                u    = Cc*xc + Dc*y
%
%  Here, dx = x - xtrim.  The user will supply dx, the initial perturbation
%  around the trim points. The output equation's J and H should also be
%  provided in Matrix form evaluated at dx. 
% 
% The opts has the following structure: 
%
% opts.xtrim = trim states for F/A-18 full nonlinear model
% opts.utrim = trim inputs for F/A-18 full nonlinear model
% opts.K     = SS Realization of the controllers used
% opts.Acl   = Closed-loop A Matrix
% opts.J     = Matrix evaluated at dx
% opts.H     = matrix evaluated at dx
% opts.Minv  = Matrix inv(M) , where M:= (I + J*Dc)
%
% Abhijit 02/16/10
% Abhijit 03/04/10

%==========================================================================
% Aircraft Physical Paramters

F18_AeroData = 1; 
f18_data; 


%==========================================================================
%-------- State ordering 

% --------- Constant Term 
V       = opts.xtrim(1); 
theta   = opts.xtrim(8); 

beta    =  x(1) + opts.xtrim(2);       % Sideslip Angle, rad
alpha   =  x(2) + opts.xtrim(3);       % Angle-of-attack, rad
p       =  x(3) + opts.xtrim(4);       % Roll rate, rad/s
q       =  x(4) + opts.xtrim(5);       % Pitch rate, rad/s
r       =  x(5) + opts.xtrim(6);       % Yaw rate, rad/s
phi     =  x(6) + opts.xtrim(7);       % Roll Angle, rad
xc      =  x(7); 

xstate  = [beta ; alpha; p; q; r; phi]; 

%--------- Input Terms
       
d_STAB   = opts.utrim(3);       % Stabilator Deflection, rad
d_RUD    = opts.utrim(2);       % Rudder Deflection, rad
d_AIL    = opts.utrim(1);       % Aileron Deflection, rad        
T        = opts.utrim(4);              % Thrust, lbf

utrim    = [d_AIL; d_RUD; d_STAB]; 

%---------- Trigonometric term

cosbeta = cos(beta); 
cos2beta3 = cos(2*beta/3);
sinbeta = sin(beta); 
tanbeta = tan(beta);  
secbeta = sec(beta); 

cosalpha = cos(alpha); 
sinalpha = sin(alpha); 

cosphi = cos(phi); 
sinphi = sin(phi); 

costheta = cos(theta); 
sintheta = sin(theta); 
sectheta = sec(theta);
tantheta =  tan(theta); 
          



%==========================================================================
% Aerodynamic Model 


% load F18AeroData; 

% % -------------------------------------------------
% Rolling Moment   
Clb     =  F18Aero.Clb_0 + F18Aero.Clb_1*alpha + F18Aero.Clb_2*alpha^2 ...
          + F18Aero.Clb_3*alpha^3  + F18Aero.Clb_4*alpha^4;
Cldr    =  F18Aero.Cldr_0 + F18Aero.Cldr_1*alpha + F18Aero.Cldr_2*alpha^2 ...
          + F18Aero.Cldr_3*alpha^3; 
Clda    =  F18Aero.Clda_0 + F18Aero.Clda_1*alpha + F18Aero.Clda_2*alpha^2 ...
          + F18Aero.Clda_3*alpha^3 ; 
Clp     =  F18Aero.Clp_0 + F18Aero.Clp_1*alpha; 
Clr     =  F18Aero.Clr_0 + F18Aero.Clr_1*alpha + F18Aero.Clr_2*alpha^2;

% % ----- Total Rolling Moment
% C_l     =  Clb*beta + Clda* d_AIL + Cldr*d_RUD + Clp*p*b/2/V + Clr*r*b/2/V;
% 
% -------------------------------------------------
%  Yawing Moment
Cnb     = F18Aero.Cnb_0 + F18Aero.Cnb_1*alpha + F18Aero.Cnb_2*alpha^2 ; 
Cndr    = F18Aero.Cndr_0 + F18Aero.Cndr_1*alpha + F18Aero.Cndr_2*alpha^2 ...
         + F18Aero.Cndr_3*alpha^3 + F18Aero.Cndr_4*alpha^4; 
Cnda    = F18Aero.Cnda_0 + F18Aero.Cnda_1*alpha + F18Aero.Cnda_2*alpha^2 ...
         + F18Aero.Cnda_3*alpha^3 ; 
Cnr     = F18Aero.Cnr_0 + F18Aero.Cnr_1*alpha; 
Cnp     = F18Aero.Cnp_0 + F18Aero.Cnp_1*alpha; 

% ----- Total Yawing Moment
% C_n     = Cnb*beta + Cnda*d_AIL + Cndr*d_RUD + Cnr*r*b/2/V + Cnp*p*b/2/V;
    
% -------------------------------------------------
% SideForce

Cyb     = F18Aero.Cyb_0 + F18Aero.Cyb_1*alpha + F18Aero.Cyb_2*alpha^2 ; 
Cyda    = F18Aero.Cyda_0 + F18Aero.Cyda_1*alpha + F18Aero.Cyda_2*alpha^2 ...
         + F18Aero.Cyda_3*alpha^3; 
Cydr    = F18Aero.Cydr_0 + F18Aero.Cydr_1*alpha + F18Aero.Cydr_2*alpha^2 ...
         + F18Aero.Cydr_3*alpha^3; 

% -------- Total Side Force
% C_Y     = Cyb*beta + Cydr*d_RUD +  Cyda*d_AIL;
         
% -------------------------------------------------
% Pitching Moment           
Cma     =  F18Aero.Cma_0 + F18Aero.Cma_1*alpha + F18Aero.Cma_2*alpha^2; 
Cmds    =  F18Aero.Cmds_0 + F18Aero.Cmds_1*alpha + F18Aero.Cmds_2*alpha^2; 
Cmq     =  F18Aero.Cmq_0 + F18Aero.Cmq_1*alpha + F18Aero.Cmq_2*alpha^2 ...
         + F18Aero.Cmq_3*alpha^3 ;
     
% --- Total Pitching Moment
% C_m     =  Cma + Cmds* d_STAB  +  Cmq*c*q/2/V;

% -------------------------------------------------
% Lift Coefficient
CLds = F18Aero.CLds_0 + F18Aero.CLds_1*alpha+ F18Aero.CLds_2*alpha^2 ...
       + F18Aero.CLds_3*alpha^3; 
   
C_lift = (-0.0204 + 5.677*alpha - 5.4246*alpha^2 + 1.1645*alpha^3)*cos2beta3;
     %+  CLds*d_STAB;

      
% -------------------------------------------------
% Drag Coefficient
Cdds = F18Aero.Cdds_0 + F18Aero.Cdds_1*alpha+ F18Aero.Cdds_2*alpha^2 ...
       + F18Aero.Cdds_3*alpha^3;  
   
C_drag =  (-1.4994 - 0.1995*alpha + 6.3971*alpha^2 - 5.7341*alpha^3 + ....
           1.4610*alpha^4) *cosbeta + 1.5036; % + Cdds*d_STAB ;
      
% -------------------------------------------------
% Form Aerodynamic forces and moments

qbar = 1/2*rho*V^2;  % Dynamic pressure
% L =  C_l*qbar*S*b ;
% M =  C_m*qbar*S*c;
% N =  C_n*qbar*S*b;
% Y =  C_Y*qbar*S ;
% Lift = C_lift*qbar*S ;
% Drag = C_drag*qbar*S ;

% Body to Wind Axis Conversion of the Aerdynamic data
% CD_w = C_drag*cosbeta - C_Y*sinbeta;
% CY_w =  C_Y*cosbeta + C_drag*sinbeta;

%==========================================================================
% Equations of Motion
%
% References: 
% 
% (i) Determination of the stability and control derivatives of the 
%     NASA F/A-18 HARV using flight data; 
%     Marcello R. Napolitano and Joelle M. Spagnuolo,NASA CR-194838, 1993. 


%----------------------------------------
% xdot = F(x) + G(x)*u
% y    = H(x) + J(x)*u

%------------------------------ betadot
Fbetad = qbar*S*(Cyb*beta*cosbeta + C_drag*sinbeta)/m/V  + p*sinalpha - r*cosalpha...
        + g*costheta*sinphi*cosbeta/V + sinbeta*(g*cosalpha*sintheta...
        - g*sinalpha*cosphi*costheta + T*cosalpha/m)/V; 
    
Gbetad =  qbar*S*[Cyda*cosbeta  Cydr*cosbeta  Cdds*sinbeta]/m/V ; 

% betadx =  qbar*S*CY_w/m/V + p*sinalpha - r*cosalpha + ....
%          g*costheta*sinphi*cosbeta/V + sinbeta*(g*cosalpha*sintheta...
%          - g*sinalpha*cosphi*costheta + T*cosalpha/m)/V;

%------------------------------ alphadot     
Falphad = -C_lift*qbar*S*secbeta/m/V + q - tanbeta*(p*cosalpha + r*sinalpha)...
          + g*(cosphi*costheta*cosalpha + sintheta*sinalpha)*secbeta/V...
          -(T*secbeta/m/V)*sinalpha;
Galphad = [0   0    -CLds*qbar*S*secbeta/m/V]; 

% alphadx = -Lift*secbeta/m/V + q - tanbeta*(p*cosalpha + r*sinalpha)...
%           + g*(cosphi*costheta*cosalpha + sintheta*sinalpha)*secbeta/V...
%           -(T*secbeta/m/V)*sinalpha;

%------------------------------ pdot
F_Cl =  qbar*S*b*(Clb*beta + Clp*p*b/2/V + Clr*r*b/2/V);
F_Cn =  qbar*S*b*(Cnb*beta + Cnr*r*b/2/V + Cnp*p*b/2/V);
Gp_da =  qbar*S*b*(Izz*Clda + Ixz*Cnda);
Gp_dr =  qbar*S*b*(Izz*Cldr + Ixz*Cndr);

Fpd = (Izz* F_Cl + Ixz* F_Cn - (Ixz*(Iyy-Ixx-Izz)*p + ...
     (Ixz^2 + Izz*(Izz - Iyy))*r)*q)/(Ixx*Izz -Ixz^2);
 
Gpd = [Gp_da  Gp_dr  0]/(Ixx*Izz -Ixz^2);

%------------------------------qdot
F_Cm = qbar*S*c*(Cma + Cmq*c*q/2/V);
G_Cm = qbar*S*c*(Cmds); 

Fqd = (F_Cm + (Izz -Ixx)*p*r + (r^2 -p^2)*Ixz)/Iyy;
Gqd = [ 0  0  G_Cm]/Iyy;

%------------------------------rdot
% rd =  ((Ixz * L  + Ixx * N + (Ixz * (Iyy - Ixx -Izz) *r +...
%       (Ixz^2 + Ixx*(Ixx - Iyy))*p)*q)/(Ixx*Izz -Ixz^2));

Gr_da =  qbar*S*b*(Ixz*Clda + Ixx*Cnda);
Gr_dr =  qbar*S*b*(Ixz*Cldr + Ixx*Cndr);

Frd = ((Ixz * F_Cl   + Ixx * F_Cn  + (Ixz * (Iyy - Ixx -Izz) *r +...
      (Ixz^2 + Ixx*(Ixx - Iyy))*p)*q)/(Ixx*Izz -Ixz^2));
  
Grd = [Gr_da  Gr_dr  0]/(Ixx*Izz -Ixz^2);
%--------------------------------------------------------------------------
% Kinetics Equation

Fphid = p + (q*sinphi + r*cosphi)*tantheta;

Gphid = [0 0 0]; 



%--------------------------------------------------------------------------
% Group F and G terms

F = [Fbetad ; Falphad; Fpd; Fqd; Frd; Fphid]; 
G = [Gbetad ; Galphad; Gpd; Gqd; Grd; Gphid]; 

%--------------------------------------------------------------------------
% Form Output Equations

H  = opts.H; 
J  = opts.J;
Minv = opts.Minv; 
%--------------------------------------------------------------------------
% Form Controller

Ac = opts.K.A;  Bc = opts.K.B;  Cc = opts.K.C;  Dc = opts.K.D; 

% -------------------------------------------------------------------------
% Form Closed Loop Model

xcldot = [F-G*Cc*xc-G*Dc*Minv*(H-J*Cc*xc ) + G*utrim; ...
          Ac*xc+Bc*Minv*(H-J*Cc*xc  )];

xcldot2 = xcldot - opts.Acl*x; 
      
% xcldot = [F-G*Cc*xc-G*Dc*Minv*(H-J*Cc*xc) ; ...
%           Ac*xc+Bc*Minv*(H-J*Cc*xc )];
