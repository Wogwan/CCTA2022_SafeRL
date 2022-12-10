% function opt = roaoptions(f,x,Name1,Value1,Name2,Value2,...)
%
% DESCRIPTION
%   Creates an options object for ROAOPTS. This can be used to estimate
%   region-of-atttraction (ROA) of polynomial system xdot = f(x), where
%   x is state vector in pvar and f(x) is the polynomial vector field.
%
% INPUTS
%      -f: polynomial vector field for which ROA is to be estimated
%      -x: state vector as in pvar
%
%   Name, Value: The options property specified by the character string
%      Name is set to Value. The setable options properties are specified
%      below.  The allowable values for each property are specified in
%      brackets with the default choice in braces.
%
%      -zV: Vector of monomials for Lyapunov function. 
%           Specifically, zV is a column vector of monomials used to 
%           specify V(x) in the vector form, V(x)=zV(x)'*C.
%           [default = monomials(x,2)]
%      -Vin: Initial Lyapunov function [ default = [] ]
%      -p: Shape factor [ default = x'*x ]
%      -z1: Monomials for s1 multipliers in beta-s1 step of V-s iteration.
%           Specifically, z1 is a column vector of monomials used to 
%           specify s1(x) in the Gram matrix form, s1(x)=z1(x)'*C*z1(x).
%           [default =  monomials(x, 0:ceil((Vdeg-p.maxdeg)/2))]
%      -z2: Monomials for s2 multiplier in gamma-s2 step of V-s iteration
%           Specifically, z2 is a column vector of monomials used to 
%           specify s2(x) in the Gram matrix form, s2(x)=z2(x)'*C*z2(x).
%           [default =  monomials(x, 1:2) ]
%      -L1: slight positive quantity to enforce positive definiteness of V
%          [ default = 0]
%      -L2: slight positive quantity to enforce negative definiteness of
%           gradient of V   [ default = 0]
%      -Q: For estimating Lyapuov function for linearized system. see linstab.m
%      -NstepBis: Number of bisection step for V-s iteration. [Default = 1]
%      -NstepMinTrace: Number of trace minimization steps for ROA iteration. 
%                   [Default = 30] 
%      -gammamax: Maximum value of gamma (Sublevel set of Lyapunov function)
%                [Default = 100]
%      -betamax: Maximum value of beta (Sublevel set of shape factor)
%                [Default = 100]
%      -sosopt: Options for V step in V-s (bisection) iteration or other 
%               steps in TraceMinimization iteration. See SOSOPTIONS 
%                  for more detail.  [Default = sosoptions]
%      -gsosopt: Options for bisection step (gamma step and beta step) in 
%                V-s (bisection) iteration. See GSOSOPTIONS for more detail. 
%                     [Default =  gsosoptions]
%      -display: Display options for iteration results. Allowable
%                    options are 'on' and 'off' [Defulat = 'off'];
%      -See SOSOPTIONS / GSOSOPTIONS for more details
%
% OUTPUT
%   opt: roaoptions object
%
% SYNTAX
%   opt = roaoptions(f,x)
%     Creates an gsosoptions object initialized with default values.
%   opt = roaoptions(f,x,Name1,Value1,Name2,Value2,...)
%     Creates an gsosoptions object with options specified by Name set
%     to the values specified by Value.
%
% See also roaoptions, L2gainoptions, gsosoptions, sosoptions

%
% 11/17/2010 Abhijit  Initial Coding

classdef roaoptions
    
    properties
        f;
        x;
        Vin;
        p ;
        zV;
        z1;
        z2;
        L2;
        L1;
        Q ;
        NstepBis = 1;
        NstepMinTrace = 30;
        gammamax = 100;
        betamax = 100;
        gsosopts = gsosoptions;
        sosopts = sosoptions;
        display ='off';
    end
    
    methods
        % Constructor
        function opt = roaoptions(varargin)
            % Check # of input args
            nin = nargin;
            
            f = varargin{1};
            x = varargin{2};
            opt.f = f;
            opt.x = x;
            
            
            if ceil(nin/2)~=floor(nin/2)
                errstr1 = 'ROAOPTIONS must have an even number of inputs';
                errstr2 = ' with Name/Value pairs specified together.';
                error([errstr1 errstr2]);
            end
            
            
            % Set Name/Value pairs:
            % Rely on default error if Name is not a public property
            for i1= 2:(nin/2)
                Name = varargin{2*(i1-1)+1};
                Value = varargin{2*i1};
                opt.(Name) = Value;
            end
            
            % initialize zV
            if isempty(opt.zV)
                opt.zV = monomials(x,2);
            end
            Vdeg = opt.zV.maxdeg;
            
            % initialize p z1 z2 L1 L2
            if isempty(opt.p)
                opt.p= x'*x;
            end
            
            if isempty(opt.L1)
                if isempty(opt.L2)
                    % This assumes that L2 will be set to be strictly pos.
                    % def.  If L2>0 and L1=0 then any feasible V wil be >0.
                    opt.L1 = 0;
                else
                    opt.L1= 1e-6*x'*x;
                end
            end
            
            if isempty(opt.L2)
                opt.L2= 1e-6*x'*x;
            end
            
            if isempty(opt.z1)
                % XXX More intelligent selection?
                opt.z1 =  monomials(opt.x, 0:ceil((Vdeg-opt.p.maxdeg)/2));
            end
            
            if isempty(opt.z2)
                % XXX More intelligent selection?
                opt.z2 =  monomials(opt.x, 1:2);
            end
            
            if isempty(opt.Q)
                opt.Q = eye(length(opt.x));
            end
        end
        
        % Set: zV
        function opt = set.zV(opt,value)
            if ismonom(value)
                opt.zV = value;
            else
                error('zV must be a vector of monomials.');
            end
        end
        
        
        % Set: Vin
        function opt = set.Vin(opt,value)
            if isa(value,'polynomial') || isempty(value)
                opt.Vin = value;
            else
                error('Lyapunov function must be a polynomial.');
            end
        end
        
        
        % Set: p
        function opt = set.p(opt,value)
            if isa(value,'polynomial')
                opt.p = value;
            else
                error('Shape function must be a polynomial.');
            end
        end
        
        
        % Set: z1
        function opt = set.z1(opt,value)
            if (ismonom(value) || isa(value,'double'))
                opt.z1 = value;
            else
                error('Multiplier of betastep must be a monomial or constant.');
            end
        end
        
        % Set: z2
        function opt = set.z2(opt,value)
            if  ismonom(value) && ~isa(value,'double')
                opt.z2 = value;
            else
                error('Multiplier of gammastep must be a monomial and non-constant.');
            end
        end
        
        % Set: L2
        function opt = set.L2(opt,value)
            if isa(value,'polynomial') || isa(value,'double') %&& isequal(value.varname,opt.x)
                opt.L2 = value;
            else
                error('L2 of gammastep must be polynomial in states.');
            end
        end
        
        % Set: L1
        function opt = set.L1(opt,value)
            if isa(value,'polynomial') || isa(value,'double')
                opt.L1 = value;
            else
                error('L1 of betastep must be a polynomial in states .');
            end
        end
        
        % Set: Q
        function opt = set.Q(opt,value)
            if isa(value,'double') && all(size(value) == [length(opt.x) length(opt.x)])
                opt.Q = value;
            else
                error('Q matrix must be a matrix of appropiate size');
            end
        end
        
        % Set: NstepBis
        function opt = set.NstepBis(opt,value)
            if isa(value,'double') && value >= 1
                opt.NstepBis = value;
            else
                error('Number of iteration steps on bisection must be a double  scalar and greater than one.');
            end
        end
        
        % Set: NstepMinTrace
        function opt = set.NstepMinTrace(opt,value)
            if isa(value,'double') && value >= 0
                opt.NstepMinTrace = value;
            else
                error('Number of iteration steps on Trace Minimization must be a double scalar.');
            end
        end
        
        % Set: sosopts
        function opt = set.sosopts(opt,value)
            if isa(value,'sosoptions')
                opt.sosopts = value;
            else
                error('sosopts must be an sosoptions object.');
            end
        end
        
        % Set: gsosopts
        function opt = set.gsosopts(opt,value)
            if isa(value,'gsosoptions')
                opt.gsosopts = value;
            else
                error('gsosopts must be a gsosoptions object.');
            end
        end
        
        % Set: gammamax
        function opt = set.gammamax(opt,value)
            if isa(value,'double') && isscalar(value) && value>0
                opt.gammamax = value;
            else
                error('gammamax must be a positive scalar.');
            end
        end
        
        % Set: betamax
        function opt = set.betamax(opt,value)
            if isa(value,'double') && isscalar(value) && value>0
                opt.betamax = value;
            else
                error('betamax must be a positive scalar.');
            end
        end
        
        % Set: display
        % 'pcontain' is an undocumented display option used by PCONTAIN
        % to display the objective function g:=-t.
        function opt = set.display(opt,value)
            AllowableVal = {'on'; 'off'};
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.display = value;
            else
                error('display can be ''on'' or ''off''. ');
            end
        end
        
        
        
    end % methods
    
end % classdef