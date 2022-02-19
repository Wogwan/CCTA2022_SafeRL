% function opt = L2reachoptions(f,x,w,p,Name1,Value1,Name2,Value2,...)
%
% DESCRIPTION
%   Creates an options object for L2REACHEST. This is used to estimate the 
%   L2 reachability ellipsoid of the following dynamics as explained in 
%   L2reachest:
%               xdot = f(x,w);
%
% INPUTS
%      -f: Vector field of polynomial system  (Ns-by-1 polynomial) 
%      -x: State  (Ns-by-1 vector of pvars)
%      -w: Input  (Nw-by-1 vector of pvars)
%      -p: Shape factor
%
%   Name, Value: The options property specified by the character string
%      Name is set to Value. The setable options properties are specified
%      below.  The allowable values for each property are specified in
%      brackets with the default choice in braces.
%
%      -zV: Vector of monomials for storage function. 
%                [default = monomials(x,1)]
%      -Vin: Initial storage function [ default = [] ]
%      -z1: Monomoails for s1 multipliers in R-s1 step of V-s iteration
%           [default =  monomials([x;w], 0:1) ]
%      -z2: Monomoails for s2 multiplier in gamma-s2 step of V-s iteration
%           [default =  monomials([x;w], 1)] 
%      -Nstep: Number of iterations for V-s iteration. [Default = 30]
%      -reachdisplay: Display options for V-s iteration results. Allowable 
%                    options are 'on' and 'off' [Defulat = 'off'];
%      -sosopt: Options for sosopt problem in iterations. See SOSOPTIONS 
%                  for more detail.  [Default = sosoptions]
%      -gsosopt: Options for gsosopt problem in V-s (Bisection) iteration. 
%                See GSOSOPTIONS for more detail. [Default =  gsosoptions] 
%      -Rmax: Minimum objective value for R-s1 step in V-s iteration is 
%             set by Rmax; by setting "gopt.minobj = -Rmax^2" 
%             [Default = 100];
%      -L1: slight positive quantity to enforce positive definiteness of V
%          [ default = 0]
%      -See SOSOPTIONS / GSOSOPTIONS for more details
%
% SYNTAX
%   opt = L2reachoptions(f,x,w,p)
%     Creates an L2reachoptions object initialized with default values.
%   opt =  L2reachoptions(f,x,w,p,Name1,Value1,Name2,Value2,...)
%     Creates an  L2reachoptions object with options specified by Name set
%     to the values specified by Value.
%
% See also roaoptions, L2gainoptions, gsosoptions, sosoptions
%

% 11/17/2010 Abhijit  Initial Coding

classdef L2reachoptions 
    
    properties  
        f;
        x;
        w;
        p;
        zV;
        Vin;
        z1;
        z2;
        Nstep = 30; 
        reachdisplay = 'off'; 
        sosopt =sosoptions; 
        gsosopt = gsosoptions ;
        Rmax = 100; 
        L1; 

    end
    
    methods
        % Constructor
           function opt = L2reachoptions(varargin)
            % Check # of input args
            nin = nargin;
            
            f = varargin{1};
            x = varargin{2};
            w = varargin{3};
            p = varargin{4}; 
            
            opt.f = f;
            opt.x = x;
            opt.w = w;
            opt.p = p;
            
            if ceil(nin/2)~=floor(nin/2)
                errstr1 = 'L2REACHOPTIONS must have an even number of inputs';
                errstr2 = ' with Name/Value pairs specified together.';
                error([errstr1 errstr2]);
            end
             
            
            %----------------- Set Name/Value pairs: 
            % Rely on default error if Name is not a public property
            for i1= 3:(nin/2)
                Name = varargin{2*(i1-1)+1};
                Value = varargin{2*i1};
                opt.(Name) = Value;
            end
            
            %----------------- initialize
            
            % initialize z1
            if isempty(opt.z1)
                opt.z1 =  monomials([x;w], 0:1);
            end

            % initialize z2
            if isempty(opt.z2)
                opt.z2 =  monomials([x;w], 1);
            end
            
            % initialize zV
            if isempty(opt.zV)
                opt.zV = monomials(x,1);
            end
            
            if isempty(opt.L1)
                opt.L1 =  0;
            end
                        
        end

        %---------------- Set: Vin
        function opt = set.Vin(opt,value)
            if isa(value,'polynomial') 
                opt.Vin = value;
            else
                error('Lyapunov function must be a polynomial.');
            end
        end

        % Set: z1
        function opt = set.z1(opt,value)
             if (ismonom(value) || isa(value,'double'))
                opt.z1 = value;
           else
                error('Multiplier s1 of R-s must be a monomial .');
            end
        end
        
        % Set: z2
        function opt = set.z2(opt,value)
             if (ismonom(value) || isa(value,'double'))
                opt.z2 = value;
           else
                error('Multiplier s2 of R-s step must be a monomial .');
            end
        end

        
        % Set: Nstep
        function opt = set.Nstep(opt,value)
            if isa(value,'double') && value >= 1
                opt.Nstep = value;
            else
                error('Number of iteration steps on bisection must be a double  scalar and greater than one.');
            end
        end
        
        % Set: Rmax
        function opt = set.Rmax(opt,value)
            if isa(value,'double') && value > 0
                opt.Rmax = value;
            else
                error('Rmax must be a positive scalar.');
            end
        end
        
%         % Set: NstepMinTrace
%         function opt = set.NstepRmax(opt,value)
%             if isa(value,'double') && value >= 0
%                 opt.NstepMinTrace = value;
%             else
%                 error('Number of iteration steps on R Maximization must be a double scalar.');
%             end
%         end
     
        % Set: sosopts
        function opt = set.sosopt(opt,value)
            if isa(value,'sosoptions')
                opt.sosopt = value;
            else
                error('sosopts must be an sosoptions object.');
            end
        end
        
        % Set: gsosopts
        function opt = set.gsosopt(opt,value)
            if isa(value,'gsosoptions')
                opt.gsosopt = value;
            else
                error('gsosopts must be a gsosoptions object.');
            end
        end
        
        % Set: L1
        function opt = set.L1(opt,value)
            if isa(value,'polynomial') || isa(value,'double')
                opt.L1 = value;
            else
                error('L1 must be a polynomial in states .');
            end
        end
        
        % Set: reachdisplay
        % 'pcontain' is an undocumented display option used by PCONTAIN
        % to display the objective function g:=-t.
        function opt = set.reachdisplay(opt,value)
            AllowableVal = {'on'; 'off'};
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.reachdisplay = value;
            else
                error('display can be ''on'' or ''off''. ');
            end
        end
                    
    end % methods
    
end % classdef