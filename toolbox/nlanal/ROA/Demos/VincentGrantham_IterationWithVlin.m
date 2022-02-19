%------------------------------------------------------------------
% Vincent & Grantham Example with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for Vincent and
%   Grantham example using the V-S iteration.  The iteration is 
%   initialized with the Lyapunov function from linearized analysis.
%
% References: 
%  1) T. Vincent and W. Grantham. Nonlinear and Optimal Control Systems.
%     Wiley, 1997.
%  2) W. Tan and A. Packard. Stability region analysis using polynomial
%     and composite polynomial Lyapunov functions and sum of squares
%     programming. IEEE Transactions on Automatic Control. 
%     53:2:565-571, 2008.
%
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2;
x = [x1;x2];
x1dot = x2;
x2dot = -(1-x1^2)*x1-x2;
f = [x1dot; x2dot];

% Create shape function and vector of Lyapunov function monomials
% The shape function is elongated along one direction.
Vdeg = 6;
th = 0.75*pi/4; 
M = [cos(th) sin(th); -sin(th) cos(th)];
p = x'*M'*diag([1 0.001])*M*x;
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

% Plot sublevel sets
domain = [-2.5 2.5 -4 4];
    
npts = 150;
[C,ph(1)]=pcontour(V,g,domain,'r',npts); hold on;
[C,ph(2)]=pcontour(p,b,domain,'b',npts); hold off;
title([' beta = ' num2str(b)]);
axis(domain)
   
% Overlay trajectories of system and additional contours of V
simdomain=[-2.5 2.5 -4 4];
Npts = 20;
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

