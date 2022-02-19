%
% This file investigates the structure of nonlinearities in the F/A-18
% Model. All the trig functions will be approximated as Taylor series for
% this investigation purpose . This is only to see where the nonlinearity is
% dominant. 
%
% Abhijit 03/04/2010
%

%==========================================================================
% Aircraft Physical Paramters

F18_AeroData = 1; 
f18_data; 

% ------------- Cases for Trimming Nonlinear Model 
 
Level       = 0; 
CoordTurn   = 1;
UnCoordTurn = 0;
opts.poly   = 0; 

% ---------- Generate Trim Values / Open-loop linear model
GenF18LinModel; 

%----------- Generate Linear Closed-Loop Models / Load SS for Controller
GenF18LinCLPModel; 



%==========================================================================
%-------- State ordering 

pvar beta alpha p q r phi xc d_AIL d_RUD d_STAB xc

% --------- Constant Term 
V       = xtrim(1); 
theta   = xtrim(8); 
T       = utrim(4); 

x  = [beta ; alpha; p; q; r; phi]; 
xclp = [x ; xc]; 

%---------- Trigonometric term

cosbeta = 1 - (beta^2)/2; 
secbeta = 1 + (beta^2)/2; 
cos2beta3 = 1 - ((2*beta/3)^2)/2;
sinbeta = beta - beta^3/6; 
tanbeta = beta + beta^3/3; 

cosalpha = 1 - alpha^2/2; 
sinalpha = alpha - alpha^3/6; 

cosphi = 1 - phi^2/2; 
sinphi = phi - phi^3/6; 

costheta = 1 - theta^2/2; 
sintheta = theta - theta^3/6; 
sectheta = 1 + theta^2/2;
tantheta =  theta + theta^3/3; 




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

Jy = jacobian(y,u);
Hy = y-Jy*u;
% ------------------------- Approximate MBinv

MB = eye(7) + Jy * CssB.D;

m11 = MB(1,1);
m12 = MB(1,2);
m13 = MB(1,3);
m71 = MB(7,1);
m72 = MB(7,2);
m73 = MB(7,3);
%---------------------- need to apprxoimate 1/x as a taylor series or LS fit.
adata = linspace(-10*d2r,70*d2r,20); 
ydata = double(subs(m11,alpha,adata));
Ydata = 1./ydata; 

pvar c0 c1 ; 
pp = c0 + c1*alpha; 
[m11inv,cfit] = pdatafit(pp,alpha,adata,Ydata); 

Minv = polynomial(eye(7));
Minv(1,1) = m11inv;
Minv(1,2) = -m12*m11inv;
Minv(1,3) = -m13*m11inv;
Minv(7,1) = -m71*m11inv;
Minv(7,2) = -(m11*m72 - m12*m71)*m11inv;
Minv(7,3) = -(m11*m73 - m13*m71)*m11inv;


% ------------------------- Approximate MRinv
MR = eye(7) + Jy * CssR.D;
MRinv = polynomial(eye(7)); 

% tmp = det(MR); 
dettmp = -MR(1,1)*MR(7,7)+  MR(1,7)*MR(7,1);  % aq - ex
adata = linspace(-10*d2r,70*d2r,20); 
ydata = double(subs(dettmp,alpha,adata));
Ydata = 1./ydata; 

pvar c0 c1 ; 
pp = c0 + c1*alpha; % + c3*alpha*beta^3 + c4*alpha^2*beta + c5*alpha^2*beta^3;

[invdetfit,cfit] = pdatafit(pp,alpha,adata,Ydata);

MRinv(1,1) = -invdetfit*MR(7,7); 
MRinv(1,2) = -invdetfit*( MR(1,7)*MR(7,2) -  MR(1,2)*MR(7,7)); 
MRinv(1,3) = -invdetfit*( MR(1,7)*MR(7,3) -  MR(1,3)*MR(7,7));  
MRinv(1,5) = -invdetfit*( MR(1,7)*MR(7,5) -  MR(1,5)*MR(7,7)); 
MRinv(1,7) =  invdetfit*MR(1,7); 

MRinv(7,1) =  invdetfit*MR(7,1); 
MRinv(7,2) = -invdetfit*( MR(7,1)*MR(1,2) -  MR(1,1)*MR(7,2)); 
MRinv(7,3) = -invdetfit*( MR(7,1)*MR(1,3) -  MR(1,1)*MR(7,3));  
MRinv(7,5) = -invdetfit*( MR(7,1)*MR(1,5) -  MR(1,1)*MR(7,5)); 
MRinv(7,7) = -invdetfit*MR(1,1);  


% ------- Form terms from the closedloop equation

GCcxcB = G*CssB.C*xc;
HJCcxcB = Hy -Jy*CssB.C*xc; 
GDcBMinv   = G*CssB.D*Minv; 
GHB       = GDcBMinv*HJCcxcB; 

GCcxcR = G*CssR.C*xc;
HJCcxcR = Hy - Jy*CssR.C*xc; 
GDcRMinv   = G*CssR.D*MRinv; 
GHR   =  GDcRMinv*HJCcxcR; 
return
% -------- Baseline
fprintf('\n -------- Baseline beta dot \n')
betadotB = [ F(1)+ GCcxcB(1)+ GHB(1)];
[g0 , g1, h] = collect(betadotB,[beta;alpha;phi])

fprintf('\n --------Baseline alpha dot \n')
alphadotB = [ F(2)+ GCcxcB(2)+ GHB(2)] 
[g0 , g1, h] = collect(alphadotB,[beta;alpha;phi])

fprintf('\n --------Baseline p dot \n')
pdotB = [ F(3)+ GCcxcB(3)+ GHB(3)] 
[g0 , g1, h] = collect(pdotB,[beta;alpha;phi])

fprintf('\n --------Baseline q dot \n')
qdotB = [ F(4)+ GCcxcB(4)+ GHB(4)]
[g0 , g1, h] = collect(qdotB,[beta;alpha;phi])

fprintf('\n --------Baseline r dot \n')
rdotB = [ F(5)+ GCcxcB(5)+GHB(5)] 
[g0 , g1, h] = collect(rdotB,[beta;alpha;phi])

fprintf('\n --------Baseline phi dot \n')
phidotB = [ F(6)+ GCcxcB(6)+ GHB(6)]
[g0 , g1, h] = collect(phidotB,[phi])


% -------- Revised
fprintf('\n --------Revised beta dot \n')
betadotR = [ F(1)+ GCcxcR(1)+ GHR(1)] 
[g0 , g1, h] = collect(betadotR,[beta;alpha;phi])

fprintf('\n --------Revised alpha dot \n')
alphadotR = [ F(2)+ GCcxcR(2)+ GHR(2)] 
[g0 , g1, h] = collect(alphadotR,[beta;alpha;phi])

fprintf('\n --------Revised p dot \n')
pdotR = [ F(3)+ GCcxcR(3)+ GHR(3)] 
[g0 , g1, h] = collect(pdotR,[beta;alpha])

fprintf('\n --------Revised q dot \n')
qdotR = [ F(4)+ GCcxcR(4)+GHR(4)]
[g0 , g1, h] = collect(qdotR,[beta;alpha;phi])

fprintf('\n --------Revised r dot \n')
rdotR = [ F(5)+ GCcxcR(5)+ GHR(5)] 
[g0 , g1, h] = collect(rdotR,[beta;alpha])

fprintf('\n --------Revised phi dot \n')
phidotR = [ F(6)+ GCcxcR(6)+ GHR(6)]
[g0 , g1, h] = collect(phidotR,[phi])
