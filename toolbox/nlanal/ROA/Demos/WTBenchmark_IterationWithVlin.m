%------------------------------------------------------------------
% Benchmark example with VS Iteration initialized with Vlin
%   Compute an estimate of the region of attraction for benchmark
%   example given in Section 3.4 of W. Tan's Thesis. The V-S iteration
%   is initialized with the Lyapunov function from linearized analysis.
%   The ROA for this example is, by construction, a known ellipsoid.
%
% References: 
%  1) W. Tan.  Nonlinear Control Analysis and Synthesis using
%     Sum-of-Squares Programming. Ph.D. Thesis, Univ. of California,
%     Berkeley, Spring 2006.
%     (See Section 3.4 for results on the benchmark example)
%
%  2) E.J. Davison and E.M. Kurak. A computational method for determining
%      quadratic Lyapunov functions for non-linear systems.  Automatica,
%      7:627-636, 1971.
%
% 11/18/2010 Updated with roaest code.  
%------------------------------------------------------------------
clear all; close all;

% Create vector field
% The ROA is given by { x: x'*B*x<1}
nx = 3;  
x=mpvar('x',nx,1);

[UU,SS,VV] = svd(randn(nx));
L = exp(2*randn(nx,1));
B = UU*diag( L )*UU';

xdot = -x + (x'*B*x)*x;
f = xdot;

% Create shape function and vector of Lyapunov function monomials
[UU,SS,VV] = svd(randn(nx));
L = exp(2*randn(nx,1));
R = UU*diag( L )*UU';
p = x'*R*x;

Vdeg = 4;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
Nsteps = 5; 

% Compute true value of beta 
% (See Sec 3.4 of W. Tan for derivation)
lmax = max(eig(B,R));
b_true = 1/lmax;

NstepMinTrace = 10;
roadisplay = 'on';

ropt = roaoptions(f,x,'zV',zV,'p',p,'z1',z1,'z2',z2,'L2',L2,...
    'NstepMinTrace',NstepMinTrace,'display',roadisplay);

% Run the iteration code
[b,V,g,s1,s2,iter] = roaest(f,x,ropt);

fprintf(' beta = %4.3f \t beta_true = %4.3f\t',b,b_true);
fprintf('beta/beta_true = %4.3f\n',b/b_true);


