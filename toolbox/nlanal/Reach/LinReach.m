function Vlin =LinReach(f,x,w,l1,sopt)
% function V = LinReach(f,x,w,l1,sopt)
%
% DESCRIPTION 
%   This function performs a linear reachability analysis for a polynomial 
%   system, xdot=f(x,w), about the equilibrium point x = 0. The linearization
%   is computed around x=0: xdot = A*x + B*w. This function generates an
%   appropiate storage function for the liinear system by solving the
%   following conditions: 
%             C1)  V - l1 >= 0
%             C2)  -grad(V)*(A*x + B*w) + w'*w >= 0 
%
% INPUTS 
%   f: Vector field of polynomial system  (Ns-by-1 polynomial)
%   x: State  (Ns-by-1 vector of pvars)
%   w: Input  (Nw-by-1 vector of pvars)
%   l1: slight positive quantity to enforce positive definiteness of V
%       [default = 0]
%   sopt: Options for solving SOS inequalities. See SOSOPTIONS for more 
%         detail.  [Default = sosoptions]
%     
% OUTPUTS
%   V: Quadratic storage function, V.  V is returned as empty if
%      the linear analysis is infeasible.
%
% SYNTAX
%   V = LinReach(f,x,w)
%   V = LinReach(f,x,w,l1)
%   V = LinReach(f,x,w,l1,sopt)
%
% EXAMPLE
%   pvar x1 x2 w;
%   x = [x1;x2];
%   f = [-2*x1+x2+x1^2 + w; x1-3*x2]
%   V = LinReach(f,x,w)
%

%  Abhijit 12/9/2010 : initial coding

% Error Checking
if nargin < 3
    error('minimum 3 input arguments needed')
end

if nargin < 4
    l1 = 0;
    sopt = sosoptions;
elseif nargin < 5
    sopt = sosoptions;
end

% Get linear systems for f
[A,B] = plinearize(f,x,w);

% Solve SOSOPT problem
V = sosdecvar('c',x);
sosc = polyconstr;
sosc(1) = V-l1 >= 0;
Vdot = jacobian(V,x)*(A*x+B*w);
sosc(2) = w'*w - Vdot >= 0;
[info,dopt]=sosopt(sosc,[x;w],sopt);

if info.feas
    Vlin = subs(V,dopt);
else
   Vlin =[];
end
