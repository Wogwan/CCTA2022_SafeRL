%------------------------------------------------------------------
% Matallana Et. Al. Example with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for Matallana, et. al
%   example using the V-S iteration.  The iteration is 
%   initialized with the Lyapunov function from linearized analysis.
%
% References: 
%  1) L.G. Matallana, A.M. Blanco, and J.A. Bandoni. Estimation of domains 
%     of attraction in epidemiological models with constant removal rates 
%     of infected individuals.  16th Argentine Bioengineering Congress 
%     and the 5th Conference of Clinical Engineering,  Journal of Physics: 
%     Conference Series,  90:1-7, 2007.
%
% 11/18/2010 Updated with roaest code.  
%------------------------------------------------------------------
clear all; close all;

% Create vector field
r = 0.6;
lambda = 0.3;
gamma = 0.8;
A = 4;
d = 0.3;
x1bar = 4.78;
x2bar = 1.78;

pvar x1 x2
x = [x1;x2];
x1dot = -d*x1 - lambda*x1*x2 - lambda*x2bar*x1 - lambda*x1bar*x2;
x2dot = -(d+gamma)*x2 + lambda*x1*x2 + lambda*x2bar*x1 + lambda*x1bar*x2;
f = [x1dot;x2dot];

% Create shape function and vector of Lyapunov function monomials
% The shape function is elongated along one direction.
Vdeg = 4;
zV = monomials(x,2:Vdeg);
th =1.5*pi/4; 
M = [cos(th) sin(th); -sin(th) cos(th)];
p = x'*M'*diag([1 0.01])*M*x;
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
NstepMinTrace = 30;
roadisplay = 'on';
ropt = roaoptions(f,x,'zV',zV,'p',p,'z1',z1,'z2',z2,'L2',L2,...
    'NstepMinTrace',NstepMinTrace,'display',roadisplay);

% Run the iteration code
[b,V,g,s1,s2,iter] = roaest(f,x,ropt);

%Extract Results

domain = [-4 4 -4 4];
npts = 150;
[C,ph(1)]=pcontour(V,g,domain,'r',npts); hold on;
[C,ph(2)]=pcontour(p,b,domain,'b',npts); hold off;
title([ 'beta = ' num2str(b)]);
axis(domain)


fprintf('Volume of {x : V(x)<=g} = %4.3f*pi\n',pvolume(V,g)/pi);

% Overlay trajectories of system and additional contours of V
simdomain=[-4 4 -4 4];
Npts = 10;
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

