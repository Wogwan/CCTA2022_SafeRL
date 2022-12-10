%------------------------------------------------------------------
% Calculation of upper bounds of the reachable sets for a 2-state dynamics
% based on the V-s iteration algorithm. The iteration is initialized
% with the storage function from linearized analysis. The example is taken
% from the following reference: 
% 
% 1) W. Tan, A. Packard, and T. Wheeler. Local gain analysis of nonlinear 
% systems. Proc. American Control Conf., pages 92{96, Minneapolis, MN, 2006.
%
%------------------------------------------------------------------

% Dynamics
pvar x1 x2 w
x = [x1;x2];
f = [-x1+x2-x1*x2^2;-x2-x1^2*x2+w];

% Parameters for estimating upper bound of reachable sets, p(x) <= beta
p = 8*x1^2-8*x1*x2+4*x2^2;
beta = 6;

% set L2reachoptions
z1 = monomials(x,0);
z2 = monomials([x;w],1);
zV = monomials(x,1);
L2reachopt = L2reachoptions(f,x,w,p,'reachdisplay','on','L1',...
      1e-6*x'*x,'zV',zV,'z1',z1,'z2',z2);

% Run iteration code to estimate upper bound of reachable sets
[R,V,s1,s2,iter] = L2reachest(f,x,w,beta,p,L2reachopt);
