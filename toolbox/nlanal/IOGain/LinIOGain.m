function [V,gammalin]= LinIOGain(f,x,u,y,l1,sopt)
% function [V,gammalin]= LinIOGain(f,x,u,y,l1,linopts)
%
% DESCRIPTION 
%   This function performs a linear input-output gain analysis for a polynomial 
%   system, xdot=f(x,w), y = h(x) about the equilibrium point x = 0. The 
%   linearization is computed around x=0: 
%                xdot = A*x + B*w.
%                 y   = C*x + D*u
%   This function generates an appropiate storage function for the liinear 
%   system by solving the following conditions: 
%     C1)  V -l1 >= 0
%     C2)  gammalin^2*u'*u - (C*x+D*u)'*(C*x+D*u) - jacobian(V,x)*(A*x+B*u) 
%
% INPUTS 
%   f: Vector field of polynomial system  (Ns-by-1 polynomial)
%   x: State  (Ns-by-1 vector of pvars)
%   u: Input  (Nu-by-1 vector of pvars)
%   y: Output (Ny-by-1 vector of pvars)
%   l1: slight positive quantity to enforce positive definiteness of V
%       [default = 0]
%   sopt: Options for solving SOS inequalities. See SOSOPTIONS for more 
%         detail.  [Default = sosoptions]
%     
% OUTPUTS
%    V: Quadratic storage function, V.  V is returned as empty if
%      the linear analysis is infeasible.
%    gammalin: Input/output gain for linear system. (Can also be calculated  
%            by infinty norm of the linear system) 
%
% SYNTAX
%   [V,gammalin]= LinIOGain(f,x,u,y)
%   [V,gammalin]= LinIOGain(f,x,u,y,l1)
%   [V,gammalin]= LinIOGain(f,x,u,y,l1,sopt)
%
% EXAMPLE
%   pvar x1 x2 u;
%   x = [x1;x2];
%   f = [-2*x1+x2+x1^2 + u; x1-3*x2];
%   y = 5*x; 
%   [V,gammalin] = LinIOGain(f,x,u,y)
%

% XXX Abhi 11/19/2010 : need to incorporate sub-optimal  

% XXX Error Checking
% sopt = linopts.linstepopt; XXX Should this even be an option
% linopts.subopt = 'off'; 


if nargin < 4
    error('minimum 4 input arguments needed')
end

if nargin < 5
    l1 = 0;
    sopt = sosoptions;
elseif nargin < 6
    sopt = sosoptions;
end

% Get linear systems for f , y
[A,B] = plinearize(f,x,u); 
[C,D] = plinearize(y,x,u); 
    
% ---- SOS Constraint to compute linearized V
pvar copt2
V = polydecvar('c',monomials(x,2),'vec');
sosineq(1) = V -l1;
sosineq(2) = copt2*u'*u - (C*x+D*u)'*(C*x+D*u) - jacobian(V,x)*(A*x+B*u);
[info,dopt]= sosopt(sosineq,[x;u],copt2,sopt);
isfeasp = info.feas;

if isfeasp
    V = subs(V,dopt);
    gammalin = sqrt(double(subs(copt2,dopt))); 
else
    V = []; gammalin =[];  
end 


