%------------------------------------------------------------------
% Pendubot example with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for the Pendubot
%   example using the V-S iteration.  The iteration is initialized
%   with the Lyapunov function from linearized analysis.
%
% References:
%  1) U. Topcu.  Quantitative Local Analysis of Nonlinear Systems.
%     Ph.D. Thesis, Univ. of California, Berkeley, Fall 2008.
%     (See Section 3.3.4 for results on this example)
%  2) U. Topcu, A. Packard, and P. Seiler, Local stability analysis
%     using simulations and sum-of-squares programming," To appear in
%     Automatica.
%
% 11/19 Erroring out in NEW and OLD Code. But, OLD code makes progress
% upto i = 26. But, new code bails out soon.
%------------------------------------------------------------------
%clear all; close all;

%error('Example in progress')

% Create vector field
pvar x1 x2 x3 x4;
x=[x1;x2;x3;x4];
x1dot = x2;
x2dot = 781.9824*x1 + 135.3145*x2 + 689.4568*x3 + 90.2337*x4;
x3dot = x4;
x4dot = 278.8839*x1*x3^2 - 1424.511*x1 - 256.587*x2 + 272.5842*x3^3 ...
    - 1249.0615*x3 - 171.1036*x4;
f = [x1dot; x2dot; x3dot; x4dot];


% Create shape function and vector of Lyapunov function monomials
M = [80.9514   14.4653   64.2197    8.3850; ...
    14.4653    3.2059   11.5393    1.8657; ...
    64.2197   11.5393   53.1873    6.7572; ...
    8.3850    1.8657    6.7572    1.0988];
p = x'*M*x;

% Run MinTrace Algorithm with Quadratic V
roadisplay = 'on';
zV = monomials(x,2);
z1 = monomials(x, 0:1 );
z2 = monomials(x, 1 );
NstepMinTrace = 50;

ropt = roaoptions(f,x,'zV',zV,'p',p,'z1',z1,'z2',z2,...
    'NstepMinTrace',NstepMinTrace,'display',roadisplay);

tic
[b,V,g,info] = roaest(f,x,ropt);
toc

% Run V-s iteration with Quartic V
% Initialize iteration using quadratic V from MinTrace iteration
Vdeg = 4;
ropt.zV = monomials(x,2:Vdeg);
ropt.z1 = monomials(x, 0:2 );
ropt.z2 = monomials(x, 1:2 );
ropt.NstepBis = 10;
ropt.Vin = V;

tic
[b,V,g,info] = roaest(f,x,ropt);
toc


fprintf('Volume of {x : V(x)<=g} = %4.3f*(4/3*pi)\n',pvolume(V,g)/(4/3*pi));

% Plot x1/x3 slice of V/p contours for ROA estimate
figure(1);
domain = [-1.5 1.5 -1.5 1.5];
V13 = subs(V,[x2;x4],[0;0]);
p13 = subs(p,[x2;x4],[0;0]);
[C,ph(1)]=pcontour(V13,g,domain,'r'); hold on;
[C,ph(2)]=pcontour(p13,b,domain,'b'); hold off;
axis(domain);

