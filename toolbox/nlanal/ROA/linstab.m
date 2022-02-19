function [V,A,P]=linstab(f,x,Q)
% function [V,A,P]=linstab(f,x,Q)
%
% DESCRIPTION 
%   This function performs a linear stability analysis for a polynomial 
%   system, xdot=f(x), about the equilibrium point x=0. The linearization
%   is computed around x=0: xdot = A*x. If the linearization is stable
%   then the Lyapunov equation is solved: A'*P+P*A = -Q.   The quadratic
%   Lyapunov function, V(x) = x'*P*x, is returned  If the linearization 
%   is not stable then V and P are returned as empty.
%
% INPUTS 
%   f: Vector field of polynomial system  (Ns-by-1 polynomial)
%   x: State  (Ns-by-1 vector of pvars)
%   Q: Postive definite matrix in the Lyapunov equation (Default = I)
%     
% OUTPUTS
%   V: Quadratic Lyapunov function, V=x'*P*x.  V is returned as empty if
%      the linearization is not stable.
%   A: Linearization of the polynomial system around x=0.
%   P: Positive definite solution of the Lyapunov equation.  This matrix
%      is used to form the quadtratic Lyapunov function. P is returned as 
%      empty if the linearization is not stable.
%
% SYNTAX
%   [V,A,P]=linstab(f,x)
%   [V,A,P]=linstab(f,x,Q)
%
% EXAMPLE
%   pvar x1 x2;
%   x = [x1;x2];
%   f = [-2*x1+x2+x1^2; x1-3*x2]
%   [V,A,P]=linstab(f,x)
%

% PJS 4/30/2009   Initial Coding

% XXX Error Checking
if nargin==2
    Nx = length(x);
    Q=eye(Nx);
end

% Linearize: xdot = A*x
A=plinearize(f,x);

% If A is stable then solve the Lyapunov equation: A'*P+P*A = -Q
ev = eig(A);
if min(real(ev)) < 0
    I = eye(size(A));
    if exist('lyap','file')==2
        P=lyap(A',Q);
    else
        P=LOCALlyap(A',Q);
    end
    V=x'*P*x;
else
    P=[];
    V=[];
end


%--------------------------------------------------------------
% Local function (Removes dependence on control toolbox lyap)
function X = LOCALlyap(A,Q)
% solves A*X + X*A' + Q = 0, using a simple
% vectorization of the Lyapunov operator, and then
% a simple Ax=b to solve for X.
 
N = size(A,1);
Avec = kron(A,eye(N)) + kron(eye(N),A')';
qvec = Q(:);
X = reshape(-Avec\qvec,[N N]);

