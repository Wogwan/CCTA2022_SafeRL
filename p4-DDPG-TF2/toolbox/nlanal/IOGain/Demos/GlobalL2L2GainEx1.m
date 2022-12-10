%--------------------------------------------------------------
% Global L2-L2 Gain Example
%--------------------------------------------------------------

%--------------------------------------------------------------
% Construct System and compute Hinf Norm of Linearization
%--------------------------------------------------------------

pvar x1 x2 u1 u2;

% Cubic nonlinear system:  xdot = f(x,u), y=C*x
u = [u1;u2];
x = [x1;x2];

Z1 = [x1;x2];
Z2 = [x1^2; x1*x2; x2^2];
Z3 = [x1^3; x1^2*x2; x1*x2^2; x2^3];
A1 = [-4 5; -1 -2];
A2 = [3 6 3; 1 2 1]/4;
A3 = [-1 0 -9 6; 0 -3 6 -7]/8;

B = [10 2; 0 1];
C = [1 1; 2 3];

xdot = (A1*Z1+A2*Z2+A3*Z3) + B*u;
y = C*x;

% Gain of Linearization
G = ss(A1,B,C,zeros(2));
[gamLin,wcFreq]=norm(G,inf)
%figure(1); sigma(G);

%--------------------------------------------------------------
% Construct an input which approximately achieves the max gain
%--------------------------------------------------------------

% The worst-case input can be constructed from the SVD of the
% system transfer function evaluated at the peak frequency
[Uw,Sw,Vw]=svd(freqresp(G,wcFreq));
wcAmp = abs(Vw(:,1));
wcPhase = angle(Vw(:,1));
twc = linspace(0,100,1e3);
twc = twc(:);
u1 = wcAmp(1)*sin(wcFreq*twc+wcPhase(1));
u2 = wcAmp(2)*sin(wcFreq*twc+wcPhase(2));
Uwc = [u1(:) u2(:)];
X0 = [0;0];
[Ywc,T,X] = lsim(G,Uwc,twc,X0);

u_L2norm = sqrt(trapz(twc,sum(Uwc.*Uwc,2)));
y_L2norm = sqrt(trapz(twc,sum(Ywc.*Ywc,2)));

% This gain should be approximately equal to gamLin
% The longer the simulation the better the approximation
gamLin_approx = y_L2norm/u_L2norm


%--------------------------------------------------------------
% Use sosopt to compute an upper bound on the gain 
% of the NL system
%--------------------------------------------------------------
sosconstr = cell(3,1);

% Define gsqinv:=1/gamma^2 as a decision variable
pvar gsqinv;

% Define variable for quadratic storage function
V = polydecvar('c',[x1; x2],'mat');

% Constraint 1 : V(x) >= 0
sosconstr{1} = V;

% Constraint 2 : gsqinv >= 0
sosconstr{2} = gsqinv;

% Constraint 3: u^Tu - y^Ty/gamma^2 - Vdot >=0
Vdot = jacobian(V,x)*xdot;
sosconstr{3} = u'*u - gsqinv*y'*y - Vdot;

% Set Objective: Min -gsqinv = Minimize gamma
obj = -gsqinv;

% Solve with sosopt
[info,dopt,sossol] = sosopt(sosconstr,[x;u],obj);
Vsol = subs(V,dopt)
gamNL = subs(gsqinv,dopt);
gamNL = double(gamNL)^(-1/2)

% Verify Conditions
%findsos(Vsol);
%Vsoldot = jacobian(Vsol,x)*xdot;
%findsos(u'*u - gsqinv*y'*y - Vsoldot);

%--------------------------------------------------------------
% Compute lower bound on gamNL using the worst-case input
% for the linear system
%--------------------------------------------------------------

% Simulate simulink model of NL system with Uwc
%Uwc = gg*Uwc;  % Check effect of scaling Uwc on NL gain
sim('GlobalL2L2GainModel1',twc(end));

uNL_L2norm = sqrt(trapz(t,sum(Unl.*Unl,2)));
yNL_L2norm = sqrt(trapz(t,sum(Ynl.*Ynl,2)));

% This gain should be approximately equal to gamLin
% The longer the simulation the better the approximation
gamNL_lower = yNL_L2norm/uNL_L2norm
