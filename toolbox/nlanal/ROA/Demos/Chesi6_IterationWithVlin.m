%------------------------------------------------------------------
% Chesi example 6 with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for Chesi 
%   example 6 using the V-S iteration.  The iteration is initialized
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
clear all; close all;

% Create vector field
pvar x1 x2 x3;
x = [x1;x2;x3];
x1dot = -x1+x3+x1^3+x1^2*x2+x1^2*x3;
x2dot = -2*x2+x3;
x3dot = x1-2*x3;
f = [x1dot; x2dot;x3dot];


% Create options for estimating region of attraction
Vdeg = 2;
p = x'*x;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
NstepMinTrace = 20;
roadisplay = 'on';

ropt = roaoptions(f,x,'zV',zV,'p',p,'z1',z1,'z2',z2,'L2',L2,...
    'NstepMinTrace',NstepMinTrace,'display',roadisplay);

% Run the iteration code
[b,V,g,s1,s2,iter] = roaest(f,x,ropt);

% extract results
fprintf('Volume of {x : V(x)<=g} = %4.3f*(4/3*pi)\n',pvolume(V,g)/(4/3*pi));

% Plot 3d V/p contours for ROA estimate
figure(1)
domain = [-10 10 -10 10 -10 10];
npts = 50*[1 1 1];
ph1 = patch(pcontour3(V,g,domain,npts));
set(ph1, 'FaceColor', 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'red' );
ph2 = patch(pcontour3(p,b,domain,npts));
xlabel('x1');ylabel('x2');zlabel('x3');
set(ph2, 'FaceColor', 'blue', 'EdgeColor', 'none' );
view(90,0)  % Set view for x2/x3 slice
ax = axis;

% Plot 2d slice of V/p contours for ROA estimate
figure(2);
V23 = subs(V,x1,0);
p23 = subs(p,x1,0);
[C,ph(1)]=pcontour(V23,g,domain(3:end),'r'); hold on;
[C,ph(2)]=pcontour(p23,b,domain(3:end),'b'); hold off;
axis(ax(3:end));

