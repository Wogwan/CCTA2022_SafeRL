%--------------------------------------------------------------
% Global Lyapunov Stability Example
%  Use SOS methods to compute a Lyapunov function that proves
%  global asymptotic stability of the equilibrium point x=0.
%
% References: 
%  1) P. Parrilo.  Structured Semidefinite Programs and 
%     Semialgebraic Geometry Methods in Robustness and Optimization.
%     Ph.D. Thesis, Cal. Institute of Technology, Spring 2000.
%     (See Section 7.2 for results on this example)
%
% 11/19/2010 Updated 
%--------------------------------------------------------------
clear all; close all;

pvar x1 x2;
x = [x1;x2];

x1dot = -x1 - 2*x2^2;
x2dot = -x2 - x1*x2 - 2*x2^3;
f = [x1dot;x2dot];

%--------------------------------------------------------------
% Use sosopt to find a Lyapunov function which proves
% x=0 is a Globally Exponentially Stable eq. pt.
%--------------------------------------------------------------

% Define decision variable for quadratic Lyapunov function
zV = monomials(x,1);
S = polydecvar('c',zV,'mat');
L1 = 1e-6*(x1^2+x2^2);
V = S+L1;

% Constraint 1 : S \in SOS
sosconstr(1) = S;

% Constraint 2: -Vdot - L2 \in SOS
% (PJS: For quadratic V, one can show that it is not possible to force
%       -Vdot-L2 \in SOS and we must settle for -Vdot \in SOS]
Vdot = jacobian(V,x)*f;
%L2 = 1e-6*(x1^2+x2^2);
%sosconstr{2} =-Vdot-L2;
sosconstr(2) =-Vdot;

% Solve with feasibility problem
[info,dopt,sossol] = sosopt(sosconstr,x);
Vsol = subs(V,dopt)

%--------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------

tmp = -4:0.5:4;
lt = length(tmp);
x0 = [tmp(:) repmat(tmp(1),[lt 1]); tmp(:) repmat(tmp(end),[lt 1]); ...
    repmat(tmp(1),[lt 1]) tmp(:); repmat(tmp(end),[lt 1]) tmp(:)];
tfinal = 10;
[xtraj,xconv]=psim(f,x,x0',tfinal);

for i1=1:length(xtraj)
    plot(xtraj{i1,2}(:,1),xtraj{i1,2}(:,2),'b'); hold on;
end
ax = [min(tmp) max(tmp) min(tmp) max(tmp)];
[C,h]=pcontour(Vsol,[0.1 0.5 1 2 5],ax,'r');
clabel(C,h,'FontSize',15,'color','k');
axis(ax);
