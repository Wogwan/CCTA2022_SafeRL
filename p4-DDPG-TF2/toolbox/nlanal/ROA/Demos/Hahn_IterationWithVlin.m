%------------------------------------------------------------------
% Hahn's example with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for Hahn's
%   example using the V-S iteration.  The iteration is initialized
%   with the Lyapunov function from linearized analysis.
%
% References: 
%  1) W. Hahn. Theory and application of Liapunov's direct method.
%     Prentice Hall, 1963.
%
%  2) W. Tan.  Nonlinear Control Analysis and Synthesis using
%     Sum-of-Squares Programming. Ph.D. Thesis, Univ. of California,
%     Berkeley, Spring 2006.
%     (See Section 3.1.4.2 for results on Hahn's example)
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2;
x = [x1;x2];
x1dot = -x1+2*x1^2*x2;
x2dot = -x2;
f = [x1dot; x2dot];

% Plot exact stability region
x1true = linspace(0.1,6);
figure(1); hold off;
plot(x1true,1./x1true,'g',-x1true,-1./x1true,'g');
axis([-6 6 -6 6]);
xlabel('x1');
ylabel('x2');
hold on;


% Create options for estimating region of attraction
Vdeg = 6;
zV = monomials(x, 2:Vdeg); 
p = x'*[14.47 18.55; 18.55 26.53]*x;
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
NstepMinTrace = 50;
roadisplay = 'on';
ropt = roaoptions(f,x,'zV',zV,'z1',z1,'z2',z2,'NstepMinTrace',NstepMinTrace,...
    'display',roadisplay,'L2',L2,'p',p);

% Run the iteration code
[b,V,g,s1,s2,iter] = roaest(f,x,ropt);

%Extract Results

domain = [-6 6 -6 6];
[C,ph(1)]=pcontour(V,g,domain,'r');
[C,ph(2)]=pcontour(p,b,domain,'b');
title(['beta = ' num2str(b)]);

