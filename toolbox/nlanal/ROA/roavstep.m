function [V,c] = roavstep(f,p,x,z,beta,gamma,s1,s2,L1,L2,opts)
% function [V,c] = roavstep(f,p,x,z,beta,gamma,s1,s2,L1,L2,opts)
%
% DESCRIPTION
%   This function solves the V step of the ROA iteration.  Specifically,
%   the ROA estimation problem can be formulated as:
%      max beta
%      subject to:
%            V-L1  in SOS
%           -( (V-gamma)+(beta-p)*s1 ) in SOS
%           -( (grad(V)*f+L2)+(gamma-V)*s2 ) in SOS
%
%   This function solves for a Lyapunov function V which satisfies the
%   three SOS constraints while holding all other variables fixed.
%
% INPUTS
%   f: Vector field of the system dynamics (N-by-1 polynomial vector)
%   p: Shape function (scalar polynomial)
%   x: State vector (N-by-1 vector of pvars)
%   z: Nz-by-1 column vector of monomials used to specify the Lyapunov
%         function decision variable V(x) in the vector form,
%         V(x) = c'*z(x).
%   beta: Level set of shape function (scalar double)
%   gamma: Level set of Lyapunov function (scalar double)
%   s1,s2: Set containment multipliers (scalar polynomials)
%   L1: Function used to enforce strict positive definiteness of
%          V
%   L2: Function used to enforce strict negative definiteness of
%          dV/dt
%   opts: SOS options structure.  See SOSOPTIONS.
%
% OUTPUTS
%   V: Lyapunov function satisfying the three SOS constraints.  V will
%     be returned as empty if no feasible solution is found.
%   c: Nz-by-1 vector of coefficients such that V(x) = c'*z(x)
%
% SYNTAX
%   [V,c] = roavstep(f,p,x,z,beta,gamma,s1,s2,opts)

% PJS 4/28/2009   Initial coding based on getNewVfeas


%------------------------------------------------------------------
% Error Checking
%------------------------------------------------------------------


%------------------------------------------------------------------
% Get options or set defaults for unspecified options
%------------------------------------------------------------------
if nargin==8
    L1 = []; L2=[]; opts =[];
elseif nargin==9
    L2 = []; opts = [];
elseif nargin==10
    opts = [];
elseif nargin~=11
    error('Invalid syntax: roavstep must have 10 or 11 input arguments');
end
if isempty(L1)
    if isempty(L2)
        L1 = 0;
    else
        L1= 1e-6*x'*x;
    end
end
if isempty(L2)
    L2= 1e-6*x'*x;
end
if isempty(opts)
    opts = sosoptions;
end

%------------------------------------------------------------------
% Set up and solve SOS feasibility problem
%------------------------------------------------------------------

% Create Lyapunov decision variable
[V,c] = polydecvar('c',z);

% Constraint 1: V-L1  in SOS
sosconstr = polyconstr;
sosconstr(1) = V >= L1;

% Constraint 2: -( (V-gamma)+(beta-p)*s1 ) in SOS
sosconstr(2) = (V-gamma) <= s1*(p-beta);

% Constraint 3: -( (grad(V)*f+L2)+(gamma-V)*s2 ) in SOS
Vdot = jacobian(V,x)*f;
sosconstr(3) = Vdot <= -L2 + s2*(V-gamma);

% Solve feasibility problem
[info,dopt] = sosopt(sosconstr,x,opts);

% Create output
if info.feas
    V = subs(V,dopt);
    c = subs(c,dopt);
else
    V=[];
    c=[];
end

