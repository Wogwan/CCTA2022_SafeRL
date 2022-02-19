%------------------------------------------------------------------
% Coutinho example 3 with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for Coutinho
%   example 3 using the V-S iteration.  The iteration is initialized
%   with the Lyapunov function from linearized analysis.
%
% References: 
%  1) D.F. Coutinho, C.E. de Souza, and A. Trofino, Stability 
%  Analysis of Implicit Polynomial Systems, vol 54, no 5, IEEE TAC, 
%  pp. 1012-1018, May 2009.
%
% 11/18/2010 Updated with roaest code.  
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2;
x = [x1;x2];
f = [(x2^2-1)*x1 ; (x1^2-1)*x2 + (1-x2^2)*x1];

% Create options for estimating region of attraction
Vdeg = 2;
M = randn(2,2); M = M'*M;
M = eye(2);
p = x'*M*x;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
NstepMinTrace = 10;
roadisplay = 'on';
ropt = roaoptions(f,x,'zV',zV,'z1',z1,'z2',z2,'L2',L2, ...
    'NstepMinTrace',NstepMinTrace,'display',roadisplay);

% Run the iteration code
[b,V,g,s1,s2,iter] = roaest(f,x,ropt);

% Extract Results
% plot sublevel sets
domain = [-1.5 1.5 -1.5 1.5];
[C,ph(1)]=pcontour(V,g,domain,'r'); hold on;
[C,ph(2)]=pcontour(p,b,domain,'b'); hold off;
title([' beta = ' num2str(b)]);
axis(domain)
hold on
plot(1,1,'o',1,-1,'o',-1,1,'o',-1,-1,'o')
hold off;

fprintf('Volume of {x : V(x)<=g} = %4.3f*pi\n',pvolume(V,g)/pi);

