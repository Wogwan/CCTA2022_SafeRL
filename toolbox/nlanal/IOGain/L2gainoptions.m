% function opt = L2gainoptions(f,x,w,h,Name1,Value1,Name2,Value2,...)
%
% DESCRIPTION
%   Creates an options object for L2GAINOPTS. This is used to estimate the 
%   L2 gain from input w to output y of the dynamics 
%           xdot = f(x,w) ; y = h(x).
%
% INPUTS
%      -f: polynomial vector field for which L2 Gain is to be estimated 
%      -x: state vector as in pvar
%      -w: input vector of interest
%      -h: output vector field. 
%
% Name, Value: The options property specified by the character string
%      Name is set to Value. The setable options properties are specified
%      below.  The allowable values for each property are specified in
%      brackets with the default choice in braces.
%
%      -zV: Vector of monomials for storage function. 
%                [default = monomials(x,1)]
%      -Vin: Initial Lyapunov function [ default = [] ]
%      -z1: Monomoails for s1 multipliers in iteration steps.
%           [default =  monomials([x;w], 0:1)]
%      -L1: slight positive quantity to enforce positive definiteness of V
%          [ default = 0]
%      -NstepBis: Number of iterations for V-s (Bisection)iteration. 
%                [Default = 1] 
%      -NstepRmax: Number of iterations for Rmax iteration. [Default = 30] 
%      -sosopt: Options for sosopt problem in iteration . See SOSOPTIONS 
%                  for more detail.  [Default = sosoptions]
%      -gsosopt: Options for gsosopt problem in V-s (Bisection) iteration. 
%                See GSOSOPTIONS for more detail. 
%                [Default = gsosoptions('minobj',0)]                      
%      -display: Display options for V-s iteration results. Allowable 
%                    options are 'on' and 'off' [Defulat = 'off'];
%
%      -See SOSOPTIONS / GSOSOPTIONS for more details
%
% OUTPUT
%    opt: L2gainoptions object
%
% SYNTAX
%   opt = L2gainoptions(f,x,w,h)
%     Creates an gsosoptions object initialized with default values.
%   opt = L2gainoptions(L2gainoptions,Name1,Value1,Name2,Value2,...)
%     Creates an gsosoptions object with options specified by Name set
%     to the values specified by Value.
%
% See also roaoptions, L2reachoptions, gsosoptions, sosoptions
%

% 11/17/2010 Abhijit  Initial Coding

classdef L2gainoptions 
    
    properties  
        f;
        x;
        w;
        h;
        zV;
        Vin;
        z1;
        L1 = 0; 
        NstepBis = 1;
        NstepRmax = 30; 
        sosopts = sosoptions;
        gsosopts = gsosoptions('minobj',0);
        display ='off'; 
    end
    
    methods
        % Constructor
        function opt = L2gainoptions(varargin)
            % Check # of input args
            nin = nargin;
            
            f = varargin{1};
            x = varargin{2};
            w = varargin{3};
            h = varargin{4}; 
            opt.f = f;
            opt.x = x;
            opt.w = w;
            opt.h = h;
            
            if ceil(nin/2)~=floor(nin/2)
                errstr1 = 'L2GAINOPTIONS must have an even number of inputs';
                errstr2 = ' with Name/Value pairs specified together.';
                error([errstr1 errstr2]);
            end
             
            
            % Set Name/Value pairs: 
            % Rely on default error if Name is not a public property
            for i1= 3:(nin/2)
                Name = varargin{2*(i1-1)+1};
                Value = varargin{2*i1};
                opt.(Name) = Value;
            end
            
            % initialize z1  zV
            
            if isempty(opt.z1)
                opt.z1 =  monomials([x;w], 0:1);
            end

            % initialize zV
            if isempty(opt.zV)
                opt.zV = monomials(x,1);
            end
                        
        end

        % Set: Vin
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
                error('Multiplier of Rstep must be a monomial .');
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
        
        % Set: NstepRmax
        function opt = set.NstepRmax(opt,value)
            if isa(value,'double') && value >= 0
                opt.NstepRmax = value;
            else
                error('Number of iteration steps on R Maximization must be a double scalar.');
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

        
         % Set: L1
        function opt = set.L1(opt,value)
            if isa(value,'polynomial') || isa(value,'double')
                opt.L1 = value;
            else
                error('L1 of betastep must be a polynomial in states .');
            end
        end
        
          
        % Set: iodisplay
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