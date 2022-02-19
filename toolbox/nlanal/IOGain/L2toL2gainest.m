function [R,V,s1,iter] = L2toL2gainest(f,x,u,h,gamma,L2gainopts)
% function [R,V,s1,iter] = L2toL2gainest(f,x,w,h,gamma,L2gainopts)
%
% DESCRIPTION
%   This function estimates the upper bound of the L2 to L2 gain of
%   the polynomial system, 
%             xdot = f(x,w),
%              y   = h(x)
%   about the equilibrium point x = 0. In brief, the L2-gain estimation 
%   problem can be formulated as:
%      Given gamma, 
%      max R
%      subject to:
%            V-L1  in SOS
%            -[(R^2 - V )s + (grad(V).f - u'*w + gamma^-2*h'*h) ] is SOS
%
%   The above problem can be formulated as V-s iteration algorithm as
%   proposed in [1]. The V-s iteration algorithm is bilinear in  
%   optimization variable and hence involves bisection. An alternative 
%   iteration procedure, namely Rmax iteration, is also formulated. This
%   procedure does not require bisection and hence much faster than V-s
%   iteration. Both the iteration algorithm is implemented in this file. 
%
% INPUTS
%   f: Vector field of polynomial system  (Ns-by-1 polynomial)
%   x: State  (Ns-by-1 vector of pvars)
%   w: Input  (Nw-by-1 vector of pvars)
%   h: Output (Ny-by-1 vector of pvars)
%   gamma: L2-L2 Induced gain
%   L2gainhopts: Options object for estimating upper bound of L2-L2 gain.  
%                See L2GAINOPTIONS for more details.
%
% OUTPUTS
%   R: Maximum input size corresponding to the given gamma 
%   V: Storage function corresponding to the maximum R
%   s1: Multiplier of R-s1 step corresponding to the maximum R
%   iter: Struture with fields V, R, s1, gamma, and time. iter
%   contains the results on the fields for each iteration.
%
% SYNTAX
%   [R,V,s1,iter] = L2toL2gainest(f,x,u,h,gamma,L2gainopts)
%
% See also roaest, L2reachest
%
% REFERENCES:
%      1) Balas, G., Packard, A., Seiler, P., Topcu, U., 2009. Robustness 
%       analysis of nonlinear systems, NASA Workshop Slides
%       http://www.aem.umn.edu/AerospaceControl/.

%==========================================================================
% Extract information from options

zV = L2gainopts.zV;
Vin = L2gainopts.Vin;
z1 = L2gainopts.z1;
NstepBis = L2gainopts.NstepBis; 
NstepRmax = L2gainopts.NstepRmax;  
iodisp = L2gainopts.display;
sopt = L2gainopts.sosopts; 
gopt = L2gainopts.gsosopts; 
L1 =  L2gainopts.L1; 

Vdeg = zV.maxdeg;
Nstep = NstepBis + NstepRmax;

%==========================================================================
% Initialize Storage Variable
c0 = cell(Nstep,1);
iter= struct('V',c0,'R',c0,'gamma',c0,'s1',c0,'time',c0);

%==========================================================================
% Construct Lyap function from linearization

if isempty(Vin)
    [Vlin, gammalin] = LinIOGain(f,x,u,h,L1,sopt);
    fprintf('Linear System Input-Output Gain %4.6f\n',gammalin)
    fprintf('\n')
end

%==========================================================================
% Back off on gammalin and resolve

if isempty(Vin)
    
    %     if gammalin < 1e-2
    ep = 0.01;
    %     else
    %         ep = 0;
    %     end
    
    gammalin = 1.01*gammalin + ep;
    [A,B] = plinearize(f,x,u);
    [C,D] = plinearize(h,x,u);
    VV = sosdecvar('c',x);
    Vdot = jacobian(VV,x)*(A*x+B*u);
    sosc = polyconstr;
    sosc(1) = VV>=L1;
    % XXX what if there is Du?
    sosc(2) = Vdot<=gammalin^2*u'*u-x'*C'*C*x;
    [info,dopt]=sosopt(sosc,[x;u],sopt);
    
    if info.feas
        V = subs(VV,dopt);
        V = V / gammalin^2;
    else
        V = Vlin;
        V = V / gammalin^2;
    end
else
    V = Vin; % XXX do we have to scale this
end

%==========================================================================
% Run V-s iteration

fprintf('\n---------------Beginning V-s iteration\n');
biscount = 0;

for i1=1:NstepBis;
    
    tic;
    
    %==================================================================
    % V Step
    % given R2 and s1 find V such that
    %
    % V - L1 is SOS
    % -[ (R^2 - V )s + (grad(V).f - w'*w + gamma^-2*z'*z) ] is SOS
    %==================================================================
    if i1 > 1
        V = sosdecvar('cV',zV);
        Vdot = jacobian(V,x)*f;
        sosc = polyconstr;
        sosc(1) = V - L1>=0 ;
        sosc(2) =  -( Vdot - u'*u + (1/gamma^2)*h'*h + (R^2- V)*s1 ) >=0;
        [info,dopt,sossol] = sosopt(sosc,[x;u],sopt);
        
        if info.feas
            V =subs(V,dopt);
        else
            if strcmp(iodisp,'on')
                fprintf('V-step is infeasible at iteration = %d\n',i1);
            end
            break;
        end
    end
    
    %==================================================================
    %  R2-s1 Step
    %  Given V max R^2 such that
    %
    %  s is SOS
    %  -[ (R^2 - V )s + (grad(V).f - u'*u + gamma^-2*h'*h) ] is SOS
    %==================================================================
    
    Vdot = jacobian(V,x)*f;
    [R2bnds,s1]= pcontain(Vdot - u'*u + (1/gamma^2)*h'*h ,V,z1,gopt);
    
    s1tmp = s1;
    
    if ~isempty(R2bnds)
        R = sqrt(double(R2bnds(1)));
    else
        if strcmp(iodisp,'on')
            fprintf('R-s1 step is infeasible at iteration = %d\n',i1);
        end
        break;
    end
    
    
    if strcmp(iodisp,'on')
        fprintf('iteration = %d  \t R = %4.6f \n',i1,R);
    end
    iter(i1).V = V;
    iter(i1).R = R;
    iter(i1).gamma = gamma;
    iter(i1).s1 = s1;
    iter(i1).time = toc;
    biscount = biscount+1;
end

if strcmp(iodisp,'on')
    fprintf('---------------Ending V-s iteration.\n');
end


%==========================================================================
% Run R-max iterarion

if strcmp(iodisp,'on')
    fprintf('\n---------------Beginning R-max iteration.\n');
end

Rmax = 1;
Rmaxcount = 0;
Rmon = 0;
for i1=1:NstepRmax
    tic;
    
    % max_{V,R} R
    pvar R2
    V = sosdecvar('cv',zV );
    sosc = polyconstr;
    sosc(1) = V>=L1;
    Vdot = jacobian(V,x)*f;
    sosc(2) = Vdot<=gamma^2*u'*u-h'*h-s1*(gamma^2*R2-V);
    %sosc(2) = Vdot<=u'*u-1/gamma^2*h'*h-s1*(R2-V);
    
    if Rmax
        % Max R
        [info,dopt,sossol] = sosopt(sosc,[x;u],-R2,sopt);
        if info.feas
            V = subs(V,dopt);
            Vdot = jacobian(V,x)*f;
            R2 = double(subs(R2,dopt));
            R = sqrt(R2);
        else
            if strcmp(iodisp,'on')
                fprintf('Switching to feas problem on iteration %d \n',i1);
            end
            V = iter(biscount+i1-1).V;
            Vdot = jacobian(V,x)*f;
            R = iter(biscount+i1-1).R;
            R2 = R^2;
            Rmax = 0;
        end
        
    else
        % Feas problem
        if Rmon
            ROld = iter(biscount+i1-1).R;
            sosc(3) = (R2 - ROld^2) >=0;
        end
        [info,dopt,sossol] = sosopt(sosc,[x;u],sopt);
        
        if info.feas
            V = subs(V,dopt);
            Vdot = jacobian(V,x)*f;
            R2 = double(subs(R2,dopt));
            R = sqrt(R2);
        else
            if strcmp(iodisp,'on')
                fprintf('R step is infeasible on iteration %d \n',i1);
            end
            break
        end
        
        % Max R2
        pvar R2
        sosc = polyconstr;
        sosc = Vdot<=gamma^2*u'*u-h'*h-s1*(gamma^2*R2-V);
        %sosc = Vdot<=u'*u-1/gamma^2*h'*h-s1*(R2-V);
        [info,dopt,sossol] = sosopt(sosc,[x;u],-R2,sopt);
        
        if info.feas
            R2 = double(subs(R2,dopt));
            R = sqrt(R2);
            
            % Switch back to Rmax if R is decreasing
            if i1 > 1
                ROld = iter(biscount+i1-1).R;
                if R < ROld
                    Rmon = 1;
                end
            end
            
        else
            if strcmp(iodisp,'on')
                fprintf('R step max is infeasible on iteration %d \n',i1);
            end
            break
        end        
    end
    
    % Feas s
    s1 = sosdecvar('cs',z1);
    sosc = polyconstr;
    sosc(1) = s1>=0;
    sosc(2) = Vdot<=gamma^2*u'*u-h'*h-s1*(gamma^2*R2-V);
    %sosc(2) = Vdot<=u'*u-1/gamma^2*h'*h-s1*(R2-V);
    [info,dopt,sossol] = sosopt(sosc,[x;u],sopt);
    if info.feas
        s1 = subs(s1,dopt);
    else
        if strcmp(iodisp,'on')
            fprintf('s1 feas step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    
    s1tmp = s1;
    pvar ee
    s1 = sosdecvar('cs',z1);
    sosc = polyconstr;
    sosc(1) = s1>=0;
    sosc(2) = Vdot<=gamma^2*u'*u-h'*h-s1*(gamma^2*R2-V) - ee*s1tmp;
    %sosc(2) = Vdot<=u'*u-1/gamma^2*h'*h-s1*(R2-V) - ee*s1tmp;
    [info,dopt,sossol] = sosopt(sosc,[x;u],-ee,sopt);
    if info.feas
        s1 = subs(s1,dopt);
    else
        if strcmp(iodisp,'on')
            fprintf('s1 feas step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    
    % Print results and store data
    if strcmp(iodisp,'on')
        fprintf('iteration=%d \t  R = %4.3f\n',i1,R)
    end
    iter(biscount+i1).V = V;
    iter(biscount+i1).s1 = s1;
    iter(biscount+i1).R = R;
    iter(biscount+i1).time = toc;
    Rmaxcount = Rmaxcount+1;
    
end
if strcmp(iodisp,'on')
    fprintf('---------------Ending R-max iteration.\n');
end

% Set outputs
iter(biscount + Rmaxcount+1:end) = [];
[tmp , idx ] = max([iter.R]);
R = iter(idx).R;
V = iter(idx).V;
s1 = iter(idx).s1;