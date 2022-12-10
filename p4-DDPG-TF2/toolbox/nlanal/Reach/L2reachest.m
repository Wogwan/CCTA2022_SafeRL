function [R,V,s1,s2,iter] = L2reachest(f,x,u,beta,p,L2reachopts)
% function [R,V,s1,s2,iter] = L2reachest(f,x,w,beta,p,L2reachopts)
%
% DESCRIPTION
%   This function estimates the upper bound of the reachable sets of
%   the polynomial system, xdot = f(x,w), about the equilibrium point x = 0.
%   The estimation employs the so called V-s iteration. Details on the
%   proposed algorithm can be found in Reference [1]. In brief, the
%   reachability estimation problem can be formulated as:
%      Given beta, 
%      max R
%      subject to:
%            V-L1  in SOS
%           -( (V-R^2)+(beta-p)*s1 ) in SOS
%           -( grad(V)*f -w'w +(R^2-V)*s2 ) in SOS
%
% INPUTS
%   f: Vector field of polynomial system  (Ns-by-1 polynomial)
%   x: State  (Ns-by-1 vector of pvars)
%   w: Input  (Nw-by-1 vector of pvars)
%   p: Shape factor 
%   beta: Sublevel set of shape factor, p
%   L2reachopts: Options for estimating reachability ellipsoid. see 
%                L2REACHOPTIONS for more details.
%
% OUTPUTS
%   R: Maximum input size, ||w||_2 =< R,  at the end of iteration
%   V: Storage function corresponding to the maximum R
%   s1: Multiplier of R-s1 step corresponding to the maximum R
%   s2: Multiplier of V step corresponding to the maximum R
%   iter: Struture with fields V, R, s1, s2 , and time. iter
%   contains the results on the fields for each iteration.
%
% SYNTAX
%   [R,V,s1,s2,iter] = L2reachest(f,x,w,beta,p,L2reachopts)
%
% See also roaest, L2toL2gainest
%
% REFERENCES:
%      1) Balas, G., Packard, A., Seiler, P., Topcu, U., 2009. Robustness 
%       analysis of nonlinear systems, NASA Workshop Slides
%       http://www.aem.umn.edu/AerospaceControl/.

% Abhijit 12/09/10 initial coding
%==========================================================================
% Extract information from options

zV = L2reachopts.zV;
Vin = L2reachopts.Vin;
z1 = L2reachopts.z1;
z2 = L2reachopts.z2;
Nstep = L2reachopts.Nstep;
Rmax  = L2reachopts.Rmax;
reachdisp = L2reachopts.reachdisplay;
sopt = L2reachopts.sosopt;
gopt = L2reachopts.gsosopt;
gopt.maxobj = 0;
gopt.minobj = -Rmax^2;
L1 =  L2reachopts.L1;


%==========================================================================
% Initialize Storage Variable
c0 = cell(Nstep,1);
iter= struct('V',c0,'R',c0,'s1',c0,'s2',c0,'time',c0);

%==========================================================================
% Construct Lyap function from linearization

if isempty(Vin)
    V = LinReach(f,x,u,L1,sopt);
    if isempty(V)
        fprintf('Linear analysis is not feasible')
        return;
    end
else
    V = Vin;
end


%==========================================================================
% Run V-s iteration

fprintf('\n---------------Beginning V-s iteration\n');

for i1=1:Nstep;
    
    tic;
    
    %==================================================================
    % V Step
    % given R2 and s1 find V such that
    %
    %  V-L1 is SOS
    %  [(beta - p) + (V - R^2)*s1] is SOS
    %  -[ (R^2 - V )s2 + (grad(V).f - u'*u ) ] is SOS
    %==================================================================
    
    if i1 > 1
        V = sosdecvar('cV',zV);
        Vdot = jacobian(V,x)*f;
        sosc = polyconstr;
        sosc(1) = V - L1  >= 0;
        sosc(2) = -((R^2-V)*s2 + Vdot- u'*u) >= 0;
        sosc(3) = (beta-p)-(R^2-V)*s1 >= 0;
        [info,dopt,sossol] = sosopt(sosc,[x;u],sopt);
        
        if info.feas
            V =subs(V,dopt);
        else
            if strcmp(reachdisp,'on')
                fprintf('V-step is infeasible at iteration = %d\n',i1);
            end
            break;
        end
    end
    
    %==================================================================
    %  R-s Step
    %  Given V max R^2 such that
    %
    %  Vdot <= u'*u - s2*(R^2-V)
    %  (p-beta) <= s1*(V-R^2)
    %  s1, s2 is SOS
    %==================================================================
    Vdot = jacobian(V,x)*f;
    pvar c; % c:=-R^2
    s1 = sosdecvar('cs1',z1);
    s2 = sosdecvar('cs2',z2);
    sosc = polyconstr;
    sosc(1) = (c+V)*s2 - Vdot + u'*u >= 0;
    sosc(2) = (beta-p)+(c+V)*s1 >= 0;
    sosc(3) = s1 >= 0;
    sosc(4) = s2 >= 0;
    [info,dopt,sossol] = gsosopt(sosc,[x;u],c,gopt);
    if info.feas
        s1 = subs(s1,dopt);
        s2 = subs(s2,dopt);
        R2 = -double(subs(c,dopt));
        R = sqrt(R2);
    else
        if strcmp(reachdisp,'on')
            fprintf('R-s step is infeasible at iteration = %d\n',i1);
        end
        break;
    end
      
    if strcmp(reachdisp,'on')
        fprintf('iteration = %d  \t R = %4.6f \n',i1,R);
    end
    iter(i1).V = V;
    iter(i1).R = R;
    iter(i1).s2 = s2;
    iter(i1).s1 = s1;
    iter(i1).time = toc;
end

if strcmp(reachdisp,'on')
    fprintf('---------------Ending V-s iteration.\n');
end


[tmp, idx ] = max([iter.R]);
R = iter(idx).R;
V = iter(idx).V;
s1 = iter(idx).s1;
s2 = iter(idx).s2;