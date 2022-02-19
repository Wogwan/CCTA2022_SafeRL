%--------------------------------------------------------------
% Global Lyapunov Stability Example
%  Use SOS methods to compute a Lyapunov function that proves
%  global asymptotic stability of the equilibrium point x=0
%  of a cubic system.
%
% References: 
%  1) P. Seiler.  Course notes for AEM 8495, UMN, Spring 2009.
%
% 11/19/2010 Checked: OLD Code works with new sosopt
%--------------------------------------------------------------
clear all; close all;

% Cubic Model: 
%    xdot = A1*Z1(x)+ A2*Z2(x)+ A3*Z3(x);
pvar x1 x2;
x = [x1;x2];
Z1 = [x1;x2];
Z2 = [x1^2; x1*x2; x2^2];
Z3 = [x1^3; x1^2*x2; x1*x2^2; x2^3];
A1 = [-4 5; -1 -2];
A2 = [3 6 3; 1 2 1]/4;
A3 = [-1 0 -9 6; 0 -3 6 -7]/8;
xdot = A1*Z1+A2*Z2+A3*Z3;

%--------------------------------------------------------------
% Use sosopt to find a Lyapunov function which proves
% x=0 is a Globally Exponentially Stable eq. pt.
%--------------------------------------------------------------

% Define decision variable for quadratic Lyapunov function
V = polydecvar('c',[x1^2; x1*x2; x2^2],'vec');

% Constraint 1 : V(x) - L1 \in SOS
L1 = (x1^2+x2^2);
sosconstr(1) = V-L1;

% Constraint 2: -Vdot - r*V \in SOS
% r is the exponential rate of convergence
Vdot = jacobian(V,x)*xdot;
r = 3.8;
sosconstr(2) =-Vdot-r*V;

% Solve with feasibility problem
[info,dopt,sossol] = sosopt(sosconstr,x);
Vsol = subs(V,dopt)

% Verify that sos constraints are satisfied
feas1 = issos(Vsol-L1)
feas2 = issos((-jacobian(Vsol,x)*xdot-r*Vsol))

%--------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------
tmp = -2:0.5:2;
lt = length(tmp);
x0 = [tmp(:) repmat(-2,[lt 1]); tmp(:) repmat(+2,[lt 1]); ...
    repmat(-2,[lt 1]) tmp(:); repmat(+2,[lt 1]) tmp(:)];
tfinal = 10;
[xtraj,xconv]=psim(xdot,x,x0',tfinal);

for i1=1:length(xtraj)
    plot(xtraj{i1,2}(:,1),xtraj{i1,2}(:,2),'b'); hold on;
end
ax = [min(tmp) max(tmp) min(tmp) max(tmp)];
[C,h]=pcontour(Vsol,[1 5 20 40 60],ax,'r');
clabel(C,h,'FontSize',15,'color','k');
axis(ax);

