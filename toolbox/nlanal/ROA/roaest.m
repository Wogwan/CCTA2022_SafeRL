function [beta,V,gamma,s1,s2,iter] = roaest(f,x,roaopts)
% function [beta,V,gamma,s1,s2,iter] = roaest(f,x,roaopts)
%
% DESCRIPTION
%   This function estimates the lower bound of the region-of-attraction of
%   the polynomial system, xdot = f(x), about the equilibrium point x = 0.
%   In brief, the ROA estimation problem for can be formulated as:
%      max beta
%      subject to:
%            V-L1  in SOS
%           -( (V-gamma)+(beta-p)*s1 ) in SOS
%           -( (grad(V)*f+L2)+(gamma-V)*s2 ) in SOS
%
%   The above problem can be formulated as V-s iteration algorithm as
%   proposed in [1]. The V-s iteration algorithm is bilinear in  
%   optimization variable and hence involves bisection. Moreover, an alternative 
%   iteration procedure, namely MinTrace iteration, is formulated. This
%   procedure does not require bisection and hence much faster than V-s
%   iteration. Both the iteration algorithm is implemented in this file. 
%
% INPUTS
%   f: Vector field of polynomial system  (Ns-by-1 polynomial)
%   x: State  (Ns-by-1 vector of pvars)
%   roaopts: Options for estimating ROA. see ROAOPTIONS for more details.
%
% OUTPUTS
%   beta:  maximum beta achieved after the iterations
%   V: Lyapunov function corresponding to the achieved maximum beta after
%      the iterations
%   gamma: Sublevel set of Lyapunov function corresponding to the achieved 
%         maximum beta after the iterations
%   s1: Multiplier s1 corresponding to the achieved maximum beta after
%      the iterations 
%   s2: Multiplier s2 corresponding to the achieved maximum beta after
%      the iterations
%   iter: Struture with fields V, beta, gamma, s1, s2 , and time. iter
%   contains the results on the fields for each iteration.
%
% SYNTAX
%   [beta,V,gamma,s1,s2,iter] = roaest(f,x,roaopts)
%
% See also roaoptions, L2toL2gainest, L2reachest
%
% REFERENCES:
%      1) Balas, G., Packard, A., Seiler, P., Topcu, U., 2009. Robustness 
%       analysis of nonlinear systems, NASA Workshop Slides
%       http://www.aem.umn.edu/AerospaceControl/.

% 11/17/2010  Abhijit  Initial Coding


%==========================================================================
% Extract information from options
% p = [];
p  = roaopts.p ;
zV = roaopts.zV ;
z1 = roaopts.z1 ;
z2 = roaopts.z2 ;
L2 = roaopts.L2 ;
L1 = roaopts.L1 ;
Q  = roaopts.Q ;
NstepBis = roaopts.NstepBis;
NstepMinTrace = roaopts.NstepMinTrace;
sopts = roaopts.sosopts;
gopts = roaopts.gsosopts;
gopts.minobj = 0;
gammamax = roaopts.gammamax;
betamax = roaopts.betamax;
display = roaopts.display;
Vin = roaopts.Vin;
% x1 = monomials(x,1:6);

Vdeg = zV.maxdeg;
Nsteps = NstepBis + NstepMinTrace;

%==========================================================================
% Initialize Storage Variable
c0 = cell(Nsteps,1);
iter= struct('V',c0,'beta',c0,'gamma',c0,'s1',c0,'s2',c0,'time',c0);

%==========================================================================
% Run V-s iteration
fprintf('\n---------------Beginning V-s iteration\n');
biscount = 0;
for i1=1:NstepBis;
    tic;
    
    %======================================================================
    % Find V step:
    % Hold s1, s2, b, g fixed and solve the Gamma Step and Beta Step
    % for Appropiate Lyapunov function with an additional constraint
    % V-L1 in SOS
    %======================================================================
    if i1==1
        if isempty(Vin)
            % Construct Lyap function from linearization
            V=linstab(f,x,Q);
        else
            V = Vin;
        end
    else
        [V,c] = roavstep(f,p,x,zV,b,g,s1,s2,L1,L2,sopts);
        if isempty(V)
            if strcmp(display,'on')
                fprintf('V-step infeasible at iteration = %d\n',i1);
            end
            break;
        end
    end
    
    
    %======================================================================
    % Gamma Step: Solve the following problem
    % {x:V(x) <= gamma} is contained in {x:grad(V)*f <0}
    % max gamma subject to
    %                 -[grad(V)*f + (gamma - V)*s2] in SOS, s2 in SOS
    %======================================================================
    gopts.maxobj = gammamax;
    [gbnds,s2]=pcontain(jacobian(V,x)*f+L2,V,z2,gopts);
    if isempty(gbnds)
        if strcmp(display,'on')
            fprintf('gamma step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    g = gbnds(1);
    
    %======================================================================
    % Beta Step: Solve the following problem
    % {x: p(x)) <= beta} is contained in {x: V(x) <= gamma}
    % max beta subject to
    %                 -[(V - gamma) + (beta - p)*s1] in SOS, s1 in SOS
    %======================================================================
    gopts.maxobj = betamax;
    [bbnds,s1]=pcontain(V-g,p,z1,gopts);
    if isempty(bbnds)
        if strcmp(display,'on')
            fprintf('beta step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    b = bbnds(1);
    
    % Print results and store iteration data
    if strcmp(display,'on')
        fprintf('iteration = %d  \t beta = %4.6f \t gamma = %4.6f\n',i1,b,g);
    end
    iter(i1).V = V;
    iter(i1).beta = b;
    iter(i1).gamma = g;
    iter(i1).s1 = s1;
    iter(i1).s2 = s2;
    iter(i1).time = toc;
    biscount = biscount+1;
end
if strcmp(display,'on')
    fprintf('---------------Ending V-s iteration.\n');
end

%==========================================================================
% Run min trace iteration
if strcmp(display,'on')
    fprintf('\n---------------Beginning min trace iteration.\n');
end

% Compute momomials of V in Gram matrix form
zVgram = gramconstraint(0,zV');

% Run iteration
mintrace=true;
oldobj = 0;
MTcount = 0;
for i1=1:NstepMinTrace;
    tic;

    % Formulate constraints with V Gram matrix as explicit dec variable
    %if isequal(L1,0) && all(isequal(zV,monomials(x,2:Vdeg)))
%     ztmp = monomials(x,2:Vdeg);
    ztmp = monomials(x,1:Vdeg);
    if isequal(L1,0) && length(zV)==length(ztmp) && all(isequal(zV,ztmp))
        [V,csV]= sosdecvar('csv', zVgram );
        Vdot = jacobian(V,x)*f;
        s1 = sosdecvar('cs',z1);
        
        sosc = polyconstr;
        sosc(1) =  -( (Vdot+L2)+(g-V)*s2 ) >=0;
        sosc(2) =  -( (V-g)+(b-p)*s1 ) >=0 ;
        sosc(3) = s1 >= 0 ;
        sosc(4) = V >= 0 ;
    else
        V = polydecvar('cv',zV );
        Vdot = jacobian(V,x)*f;
        [sV,csV] = sosdecvar('csv', zVgram );
        s1 = sosdecvar('cs',z1);
        
        sosc = polyconstr;
        sosc(1) =  -( (Vdot+L2)+(g-V)*s2 ) >=0;
        sosc(2) =  -( (V-g)+(b-p)*s1 ) >=0 ;
        sosc(3) = s1 >= 0 ;
        sosc(4) = V-L1==sV;
        sosc(5) = sV>=0;
    end
    obj = trace(csV);
    
    % min trace of V Gram matrix
    if mintrace
        % Min trace
        [info,dopt,sossol] = sosopt(sosc,x,obj,sopts);
        if info.feas
            V = subs(V,dopt);
            Vdot = jacobian(V,x)*f;
            s1 = subs(s1,dopt);
            obj = double(subs(obj,dopt));
            
            % Check progress of objective function
            if i1>1
                tol = 0.01;
                rval= (oldobj(i1-1)-obj)/oldobj(i1-1);
                if rval<tol && mintrace==true
                    if strcmp(display,'on')
                        fprintf('Switching to feas problem on iteration %d \n',i1);
                    end
                    mintrace=false;
                end
            end
            oldobj(i1) = obj;
        else
            mintrace = false;
            if strcmp(display,'on')
                fprintf('Switching to feas problem on iteration %d \n',i1);
            end
        end
    end
    
    % Solve V/S1 as feasibility problem
    if mintrace==false
        [info,dopt,sossol] = sosopt(sosc,x,sopts);
        if info.feas
            V = subs(V,dopt);
            Vdot = jacobian(V,x)*f;
            s1 = subs(s1,dopt);
            % obj = subs(obj);
        else
            if strcmp(display,'on')
                fprintf('V/s1 feas step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        
        % Max gamma
        pvar g;
        sosc = polyconstr;
        sosc(1) =  -( (Vdot+L2)+(g-V)*s2 ) >=0;
        [info,dopt,sossol] = sosopt(sosc,x,-g,sopts);
        if info.feas
            g = double(subs(g,dopt));
        else
            if strcmp(display,'on')
                fprintf('Max gamma step infeasible at iteration = %d\n',i1);
            end
            break;
        end
    end
    
    % Max Beta
    pvar b;
    sosc = polyconstr;
    sosc(1) =  -( (V-g)+(b-p)*s1 ) >=0 ;
    [info,dopt,sossol] = sosopt(sosc,x,-b,sopts);
    if info.feas
        b = double(subs(b,dopt));
    else
        if strcmp(display,'on')
            fprintf('Max beta step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    
    % S2 Feas
    s2 = sosdecvar('cs',z2);
    sosc = polyconstr;
    sosc(1) = s2 >=0 ;
    sosc(2) =  -( (Vdot+L2)+(g-V)*s2 ) >=0;
    [info,dopt,sossol] = sosopt(sosc,x,sopts);
    if info.feas
        s2 = subs(s2,dopt);
    else
        if strcmp(display,'on')
            fprintf('s2 feas step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    
    % Print results and store data
    if strcmp(display,'on')
        fprintf('iteration = %d  \t beta = %4.6f \t gamma = %4.6f\n',i1,b,g);
    end
    iter(biscount+i1).V = V;
    iter(biscount+i1).s1 = s1;
    iter(biscount+i1).s2 = s2;
    iter(biscount+i1).beta = b;
    iter(biscount+i1).gamma = g;
    iter(i1).time = toc;    
    MTcount = MTcount+1;
end
if strcmp(display,'on')
    fprintf('---------------Ending min trace iteration.\n');
end

% Set outputs
iter(biscount + MTcount+1:end) ; 
[tmp , idx ] = max([iter.beta]);
beta = iter(idx).beta;
V = iter(idx).V;
s1 = iter(idx).s1;
s2 = iter(idx).s2;
gamma = iter(idx).gamma;


