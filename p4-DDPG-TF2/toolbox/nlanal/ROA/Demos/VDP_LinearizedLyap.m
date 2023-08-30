%------------------------------------------------------------------
% VDP with Linearized Lyapunov Function
%   Compute an estimate of the region of attraction for the van 
%   der Pol oscillator using the Lyapunov function obtained
%   via linearization.
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

% Construct Lyap function from linearization 
Q = eye(2);
%Q = diag([1 2]);
%Q = diag([5 2]);
Vlin=linstab(f,x,Q);

% Create multiplier and function to force dV/dt<0
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);

% get gamma:  maximize level set of V contained in 
%  set where dV/dt<0
opts = [];
opts = gsosoptions; 
[gbnds,s2]=pcontain(jacobian(Vlin,x)*f+L2,Vlin,z2,opts);
g = gbnds(1);
fprintf('gamma = %4.3f\n',g);
title(['gamma = ' num2str(g)])

% Plot V contour for ROA estimate 
domain = [-3 3 -3 3];
[C,ph(1)]=pcontour(Vlin,g,domain,'r');
hold off
