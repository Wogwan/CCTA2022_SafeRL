%------------------------------------------------------------------
% VDP with VS Iteration Initialized with Vball
%   Solve for V which maximizes beta subject to dV/dt<0 on p<=beta
%   Then use this V to initialize the V-S iteration
%
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2;
x = [x1;x2];
x1dot = -x2;
x2dot = x1+(x1^2-1)*x2;
f = [x1dot; x2dot];

% Plot VDP limit cycle from backward simulation
X0 = [1/2; 1/2];
Tf = 100;
[xtraj,isconv]=psim(-f,x,X0,Tf);
t = xtraj{1};
xtraj = xtraj{2};
idx = min(find(t>=Tf*0.8));

figure(1)
hold off;
plot(xtraj(idx:end,1),xtraj(idx:end,2));
xlabel('x1');
ylabel('x2');
hold on;

% Create shape function and vector of Lyapunov function monomials
Vdeg = 6;
M = randn(2); [UU,SS,VV]=svd(M);
M = UU*diag([10 1])*UU';
p = x'*M*x;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);

% Construct Lyap function by
%    max beta s.t.   -(Vdot+L2)+(p-beta) >= 0
b0 = 1.0;
V0 = polydecvar('c',zV,'vec');
s20 = polydecvar('s',z2,'vec');
sosconstr = [V0-L2; -(jacobian(V0,x)*f+L2)+(p-b0)*s20];
[info,dopt,sossol]=sosopt(sosconstr,x);
V0 = subs(V0,dopt);

% ROA Options Object
NstepMinTrace = 35;
roadisplay = 'on';
ropt = roaoptions(f,x,'zV',zV,'p',p,'z1',z1,'z2',z2,'L2',L2,...
    'NstepMinTrace',NstepMinTrace,'display',roadisplay,'Vin',V0);

% Run the iteration code
[b,V,g,s1,s2,iter] = roaest(f,x,ropt);

% plot sublevel sets of Lyapunov function and shape function
domain = [-3 3 -3 3];
[C,ph(1)]=pcontour(V,g,domain,'r');
[C,ph(2)]=pcontour(p,b,domain,'b');
title([' beta = ' num2str(b)]);


