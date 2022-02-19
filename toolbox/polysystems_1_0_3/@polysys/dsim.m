function [t,x,y] = dsim(sys,timeSpan,state0,input)
%DSIM  Simulate a discrete-time POLYSYS model.
%
%   [T,X,Y] = DSIM(SYS,TSPAN,X0,U)  simulates the discrete-time POLYSYS
%   model SYS starting from the initial condition X0.  TSPAN is a vector of
%   time points at which the solution is computed.  The input U is an
%   NT-by-NU array, where NT is the length of TSPAN and NU is the number of
%   inputs.
%
%   [T,X,Y] = DSIM(SYS,TSPAN,X0)  same as above with all inputs set to zero.
%
%   [T,X,Y] = DSIM(SYS,TSPAN)  same as above with zero initial condition.
%
%   See also POLYSYS, SIM.


error(nargchk(2,4,nargin,'struct'));
error(nargoutchk(0,3,nargout,'struct'));

if sys.sampleTime == 0
    error('POLYSYS:dsim:sysType',...
        'DSIM only works for discrete-time systems, Use SIM instead.')
end


% Check the simulation time vector.
if length(timeSpan) < 2
    error('POLYSYS:dsim:tSpanSize',...
        'Simulation time vector must have at least two points.')
end
dt = diff(timeSpan);
if not(all(dt>0))
    error('POLYSYS:dsim:tSpanDecrease',...
        'Simulation time vector must be strictly increasing.')
end
if all(dt==dt(1))
    tSample = dt(1);
else
    error('POLYSYS:dsim:tSpanSpacing',...
        'Simulation time vector must have regularly spaced entries.')
end
if sys.sampleTime > 0 && not(sys.sampleTime == tSample)
    error('POLYSYS:dsim:tSpanSampling',...
        'Simulation time vector does not match the system sampling time.')
end
timeSpan = reshape(timeSpan,length(timeSpan),1);


% Check the initial state and determine how many simulations
nStates = length(sys.states);
if nargin < 3 || isempty(state0)
    % Use default initial condition.
    state0 = zeros(nStates,1);
    nSims = 1;
elseif isnumeric(state0)
    % Determine how many simulations to perform.
    [state0Rows,state0Cols] = size( state0 );
    if state0Rows == nStates
        nSims = state0Cols;
    elseif state0Cols == nStates
        state0 = state0';
        nSims = state0Rows;
    else
        error('POLYSYS:dsim:state0Size', ...
            'Initial state array does not have proper dimensions.')
    end
else
    error('POLYSYS:dsim:state0Class', ...
        'Initial state must be a numeric vector.')
end

% Use default input if none specified.
if nargin < 4 || isempty(input)
    input = zeros(length(timeSpan),length(sys.inputs));
end

% Verify the dimensions of the system input.
if not(length(timeSpan)==size(input,1))
    error('POLYSYS:dsim:inputSamples',...
        'Input signal must be specified at each time point.')
elseif not(length(sys.inputs)==size(input,2))
    error('POLYSYS:dsim:inputDim',...
        'Input signal has the wrong number of columns.')
end
    
% Convert the system dynamics into anonymous functions.
[f,g] = local_makeFunctions(sys);

% Simulate the system.
nSteps = length(timeSpan);
nOutputs = size(sys.orMap,1);
t = cell(1,nSims);
x = cell(1,nSims);
y = cell(1,nSims);
for i = 1:nSims
    x_i = zeros(nStates,nSteps);
    y_i = zeros(nOutputs,nSteps);
    x_i(:,1) = state0(:,i);
    y_i(:,1) = g(x_i(:,1),input(1,:)');
    for j = 2:nSteps
        x_i(:,j) = f(x_i(:,j-1),input(j-1,:)');
        y_i(:,j) = g(x_i(:,j),input(j,:)');
    end
    t{i} = timeSpan;
    x{i} = x_i';
    y{i} = y_i';
end

if nSims == 1
    t = t{1};
    x = x{1};
    y = y{1};
end



function [f,g] = local_makeFunctions(sys)

% Create a function handle for the state transition map.
stMask = pevalIndex( sys.stMap, [sys.inputs;sys.states] );
stCoef = sys.stMap.coefficient;
stDeg = sys.stMap.degmat;
nStates = size(sys.stMap,1);
f = @( xIN, uIN )( reshape( ...
        peval( stMask*[uIN;xIN],stCoef,stDeg),...
    nStates,1) );

% Create a function handle for the output response map.
orMask = pevalIndex( sys.orMap, [sys.inputs;sys.states] );
orCoef = sys.orMap.coefficient;
orDeg = sys.orMap.degmat;
nOutputs = size(sys.orMap,1);
g = @( xIN, uIN )( reshape( ...
        peval( orMask*[uIN;xIN],orCoef,orDeg),...
    nOutputs,1) );