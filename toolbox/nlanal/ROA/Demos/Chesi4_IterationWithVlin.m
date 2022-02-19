%------------------------------------------------------------------
% Chesi example 4 with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for Chesi 
%   example 4 using the V-S iteration.  The iteration is initialized
%   with the Lyapunov function from linearized analysis.
%
% References: 
%  1) G. Chesi, A. Garulli, A. Tesi and A. Vicino. LMI-based
%     computation of optimal quadratic Lyapunov functions for
%     odd polynomial systems. Int. J. Robust Nonlinear Control
%     15:35–49, 2005.
%  2) U. Topcu, A. Packard, and P. Seiler, Local stability analysis 
%     using simulations and sum-of-squares programming," To appear in 
%     Automatica. 
%  3) U. Topcu.  Quantitative Local Analysis of Nonlinear Systems.
%     Ph.D. Thesis, Univ. of California, Berkeley, Fall 2008.
%     (See Section 3.3.2 for results on this example)
%  4) W. Tan and A. Packard. Stability region analysis using polynomial
%     and composite polynomial Lyapunov functions and sum of squares
%     programming. IEEE Transactions on Automatic Control. 
%     53:2:565-571, 2008.
%
%------------------------------------------------------------------
%clear all; close all;

% Create vector field
pvar x1 x2;
x = [x1;x2];
x1dot = x2;
x2dot = -2*x1-x2+x1*x2^2-x1^5+x1*x2^4+x2^5;
f = [x1dot; x2dot];

% Create options for estimating region of attraction
Vdeg = 4;
p = x'*x;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
NstepMinTrace = 1;
roadisplay = 'on';

%ropt = roaoptions(f,x,'zV',zV,'p',p,'z1',z1,'z2',z2,'L2',L2,...
%    'NstepMinTrace',NstepMinTrace,'display',roadisplay,'NstepBis',30);

ropt = roaoptions(f,x,'zV',zV,'p',p,'z1',z1,'z2',z2,'L2',L2,...
    'display',roadisplay,'NstepMinTrace',40);
ropt.L1 = 0;

% Run the iteration code
[b,V,g,s1,s2,iter] = roaest(f,x,ropt);

% plot sublevel sets of Lyapunov function and shape function
domain = [-1.5 1.5 -1.5 1.5];

[C,ph(1)]=pcontour(V,g,domain,'r'); hold on;
[C,ph(2)]=pcontour(p,b,domain,'b'); hold off;
title([ 'beta = ' num2str(b)]);
axis(domain)

fprintf('Volume of {x : V(x)<=g} = %4.3f*pi\n',pvolume(V,g)/pi);

% Overlay trajectories of system and additional contours of V
simdomain=[-0.9 0.9 -0.9 0.9];
Npts = 8;
x1ic = linspace(simdomain(1),simdomain(2),Npts);
x2ic = linspace(simdomain(3),simdomain(4),Npts);
x0 = [x1ic x1ic repmat(x1ic(1),1,Npts) repmat(x1ic(end),1,Npts);
    repmat(x2ic(1),1,Npts) repmat(x2ic(end),1,Npts) x2ic x2ic];
xtraj = psim(f,x,x0,20);

figure(1);
ax = axis;
hold on;
[C,ph(1)]=pcontour(V,g*[0.3 0.1 0.05 0.01],domain,'k',[1e3 1e3]); 
for i1=1:size(x0,2)
    plot(xtraj{i1,2}(:,1),xtraj{i1,2}(:,2),'g');hold on;
end
hold off;
axis(ax);

