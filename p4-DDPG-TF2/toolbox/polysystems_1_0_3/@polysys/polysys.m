function sys = polysys(varargin)
%POLYSYS  Dynamical system described by polynomial equations.
%
%   SYS = polysys(STMAP,ORMAP,STATES,INPUTS)  creates a continuous-
%   time system, where STMAP is a POLYNOMIAL object representing the state
%   transition map, ORMAP is a POLYNOMIAL representing the output response
%   map, STATES is a POLYNOMIAL vector of monomials indicating which
%   variables are states, and INPUTS is a POLYNOMIAL vector of monomials
%   indicating which variables are inputs.  STATES and INPUTS must be vectors 
%   of unique monomials with coefficient 1, and they cannot have any variable 
%   names in common.  STMAP and ORMAP must be column vectors of class 
%   POLYNOMIAL and STMAP must be the same size as STATES.
%
%   SYS = polysys(STMAP,ORMAP,STATES,INPUTS,TS)  creates a discrete-
%   time system with sampling time TS.  If TS == -1, the system will
%   inherit its sampling time.
%
%   SYS = polysys(OBJ)  will try to convert OBJ to a POLYSYS.  If OBJ is
%   already a POLYSYS object, then SYS = OBJ.
% 
%   SYS = polysys()  creates an empty POLYSYS object with no state
%   transition map and no output response map.
%
%   Note:
%     When a POLYSYS object is created, the arguments STATES and INPUTS
%     are used to establish the order of the state and inputs.  However, to
%     avoid naming conflicts and ambiguity, the variables are renamed so 
%     that STATES(i) is a monomial named 'xi' and INPUTS(j) is a monomial
%     named 'uj' for each appropriate i and j.  When two systems are 
%     connected (using mtimes, blkdiag, parallel, etc.), the states and 
%     inputs of the "subordinate" system are renamed so that they continue 
%     the sequence of the "superior" system's states and inputs.  For 
%     example, if SYS1.states=[x1;x2], then before connecting SYS1 to SYS2, 
%     the states of SYS2 are renamed to [x3;x4;...].
%
%   See also POLYNOMIAL.
%      interconnection methods:  polysys/blkdiag, polysys/mtimes,
%         polysys/vertcat, polysys/horzcat, polysys/plus,
%         polysys/parallel
%      simulation methods:  polysys/sim, polysys/dsim.

% 7.19.2007: TJW - Initial coding, minimal help documentation.
% 7.22.2007: TJW - Changed class name from twpds to polysys.
%                - Help documentation and notes ("See also..." needs work).
% 8.06.2007: TJW - Changed order of inputs so name is optional and last.


%% If possible, convert input to a polysys.

if nargin == 1
    
    if isa(varargin{1},'polysys')
        sys = varargin{1};
    else
        sys = toPolysys(varargin{1});
    end

%%  Create an 'empty' object with correct data types.

elseif nargin == 0
    
    sys.stMap = polynomial();
    sys.stStateIndex = [];
    sys.stInputIndex = [];
    sys.orMap = polynomial();
    sys.orStateIndex = [];
    sys.orInputIndex = [];
    sys.states = polynomial();
    sys.inputs = polynomial();
    sys.sampleTime = 0;
    sys.name = '';
    sys.isDynamic = false(1);
    sys.hasDirectFeedthrough = false(1);

    sys = class(sys,'polysys');
    superiorto('polynomial','lti','ss','tf','zpk','double');


%% Create new polysys object.

elseif nargin >= 4
    
    % Make sure first four inputs are all polynomial objects or empty arrays.
    for i = 1:4
        if not(isa(varargin{i},'polynomial'))
            if isempty(varargin{i})
                varargin{i} = polynomial();
            else
                error('POLYSYS:inputClass',...
                    'The first four inputs must be polynomial objects or an empty arrays.')
            end
        end
    end

    % Load the input data into our struct.
    sys.stMap = varargin{1};
    sys.stStateIndex = [];
    sys.stInputIndex = [];
    sys.orMap = varargin{2};
    sys.orStateIndex = [];
    sys.orInputIndex = [];
    sys.states = varargin{3};
    sys.inputs = varargin{4};
    if nargin < 5
        sys.sampleTime = 0;
    else
        sys.sampleTime = varargin{5};
    end
    if nargin < 6
        sys.name = '';
    else
        sys.name = varargin{6};
    end
    sys.isDynamic = not(isempty(sys.stMap));

    % Check the system name for consistency
    if isempty(sys.name)
        sys.name = '';
    elseif not(ischar(sys.name))  ||  size(sys.name,1) > 1
        error('POLYSYS:inputClass',...
            'System name must be a string.')
    end
  
    % Check the variables
    if not(isMonoVector(sys.states)) || not(isMonoVector(sys.inputs) )
        error('POLYSYS:varVector',...
            'States and inputs must be empty or polynomial vectors of monomials\nwith coefficient 1.');
    end
    allVars = [ sys.states; sys.inputs ];
    if length(unique(allVars.varname)) ~= length(allVars)
        error('POLYSYS:varVector',...
            'States and inputs cannot have any variables in common.')
    end
    
    % Check state-transition map (if we have one).
    if sys.isDynamic
        if size(sys.stMap,2) ~= 1
            error('POLYSYS:stMapWide',...
                'State transition map must be a column vector.')
        end
        if any(not( ismember(sys.stMap.varname, allVars.varname) ))
            error('POLYSYS:stMapExtra',...
                'State transition map contains variables that are not states or inputs.')
        end
        if length(sys.stMap) ~= length(sys.states)
            error('POLYSYS:stMapLong',...
                'State transition map must have same length as states.')
        end
    else
        if not(isempty(sys.states))
            error('POLYSYS:statesInStatic',...
                'A static map cannot have states.')
        end
    end
 
    % Check output response map.
    if isempty(sys.orMap)
        error('POLYSYS:orMapMissing',...
            'Output response map cannot be set to empty.')
    end
    if size(sys.orMap,2) ~= 1
        error('POLYSYS:orMapWide',...
            'Output response map must be a column vector.')
    end
    if any(not( ismember(sys.orMap.varname, allVars.varname) ))
        error('POLYSYS:orMapExtra',...
            'Output response map contains variables that are not states or inputs.')
    end
    
    % 'Normalize' the names of the states and inputs.
    % Need to change varnames in stMap,orMap,states,inputs
    % Need to create stStateIndex,stInputIndex,orStateIndex,orInputIndex
    if sys.isDynamic
        [sys.stStateIndex,sys.stInputIndex] = ...
            polyVarsIndex(sys.stMap,sys.states,sys.inputs);
    end
    [sys.orStateIndex,sys.orInputIndex] = ...
        polyVarsIndex(sys.orMap,sys.states,sys.inputs);
    sys = normalizeVars(sys,0,0);
        
    % See if system has direct feedthrough
    if any(ismember(sys.orMap.varname, sys.inputs.varname))
        sys.hasDirectFeedthrough = true(1);
    else
        sys.hasDirectFeedthrough = false(1);
    end
    
    % Check the sample time.
    if isempty(sys.sampleTime)
        sys.sampleTime = 0;
    else
        if isscalar(sys.sampleTime)
            if (sys.sampleTime < 0)  &&  (sys.sampleTime ~= -1)
                error('POLYSYS:badSampleTime',...
                    'Negative sample time not allowed (except -1 to mean inherited).')
            end
        else
            error('POLYSYS:badSampleTime',...
                'Sample time must be a scalar.')
        end
    end
    
    % Create the polysys object
    sys = class(sys,'polysys');
    superiorto('polynomial','lti','ss','tf','zpk','double');

    
end % End 'if nargin...'