function [t,state,output] = sim(sys,timeSpan,state0,input,options,solver)
%SIM  Simulate the dynamics of a POLYSYS object.
%
%   [T,X,Y] = sim(SYS,TSPAN,X0)  with TSPAN = [T0 TFINAL] integrates the
%   dynamics of SYS from time T0 to TFINAL with initial conditions X0 and
%   computes the output Y that corresponds to the state trajectory X.  To
%   obtain solutions at specific time T0,T1,...,TFINAL (all increasing or
%   all decreasing), use TSPAN = [T0 T1 ... TFINAL].  To do multiple
%   simulations, each starting from a different intial condition, use 
%   X0 = [X01 X02 ... X0M].  Then, T, X, and Y will be 1-by-M cell arrays,
%   where the ith cell contains the results of the simulation starting from
%   X0i.  By default, SIM sets all system inputs to zero and integrates
%   using ODE45.
%
%   [T,X,Y] = sim(SYS,TSPAN,X0,U)  simulates as above with input U.  The
%   first column of U is a vector of time points, and the remaining columns
%   are the system inputs.  Setting U = [] is the same as setting U to
%   zero.
%
%   [T,X,Y] = sim(SYS,TSPAN,X0,U,OPT,SOLVER)  uses the OPT structure
%   created by ODESET and the ODE solver specified in the string SOLVER.
%   Any of the built-in ODE solvers may be used.
%
%   See also POLYSYS, ODE45, ODESET.

% 08.??.2007: TJW - Initial coding.
% 01.19.2008: TJW - Moved the loop for multiple simulations to
%                   local_contSolver so that the odefun's only need to be
%                   created once.
% 01.20.2008: TJW - Tried to speed up by creating intermediate variables to
%                   reduce the complexity of the anonymous functions.
%                   Modified to use peval4 instead of peval, which relied
%                   heavily on reshape.
%                 - Help documentation.
% 01.22.2008: TJW - Extended sim to include discrete-time simulations.
%                   Did a quick check to see that it works but no unit tests.
% 01.24.2008: TJW - Commented out the discrete-time code and wrote dsim.m.
% 12.17.2008: TJW - Removed the usage of peval4.c

error(nargchk(2,6,nargin,'struct'));
error(nargoutchk(0,3,nargout,'struct'));

% Check the initial state.
nStates = length(sys.states);
if nargin < 3
    state0 = zeros(nStates,1);
    nSims = 1;
elseif isnumeric(state0)
    [state0Rows,state0Cols] = size( state0 );
    if state0Rows == nStates
        nSims = state0Cols;
    elseif state0Cols == nStates
        state0 = state0';
        nSims = state0Rows;
    else
        error('POLYSYS:sim:state0Size', ...
            'Initial state array does not have proper dimensions.')
    end
else
    error('POLYSYS:sim:state0Class', ...
        'Initial state must be a numeric vector.')
end

% Use default input if not specified.
if (nargin<4) 
    input = [];
end

% If system has continuous-time dynamics...
if sys.sampleTime == 0
    
    %%%% Continuous-time simulation.
    
    % Check the ODE options structure.
    if nargin < 5
        options = odeset();
    else
        options = odeset(options);
    end

    % Check the solver name.
    if nargin < 6
        solver = 'ode45';
    elseif not(ischar(solver))
        error('POLYSYS:sim:solverClass',...
            'Solver must be a string specifying the name of an ODE solver.')
    end

    % Call the continuous-time solver.
    if nargout > 2
        [t,state,output] = local_contSolver(sys,timeSpan,state0,input, ...
            options,solver,nSims);
    else
        [t,state] = local_contSolver(sys,timeSpan,state0,input, ...
            options,solver,nSims);
    end
else

    %%%% Discrete-time simulation.
    error('POLYSYS:sim:discrete',...
        'Use DSIM for discrete-time systems.')

%     % Check the interpolation option.
%     if nargin < 5
%         interpMethod = 'zoh';
%     elseif ischar(options) && ismember(lower(options),{'zoh','foh'})
%         interpMethod = options;
%     else
%         error('POLYSYS:sim:interpMethod',...
%             'Interpolation method must be either ''zoh'' or ''foh''.')
%     end
%     
%     % Call the discrete-time solver.
%     if nargout > 2
%         [t,state,output] = local_discSolver(sys,timeSpan,state0,input,...
%             interpMethod,nSims);
%     else
%         [t,state] = local_discSolver(sys,timeSpan,state0,input,interpMethod,nSims);
%     end
    
end
    
% If we only do one sim, return double arrays rather than cell arrays.
if nSims == 1
    t = t{1};
    state = state{1};
    if nargout > 2
        output = output{1};
    end
end



%% Continuous-time solver.
function [t,x,y] = local_contSolver(sys,timeSpan,initialState,input,...
    options,solver,nSims)

% Verify the input.
if isempty(input)
    input = zeros(1,length(sys.inputs)+1);
elseif isnumeric(input) 
    if size(input,2) ~= length(sys.inputs)+1
        error('POLYSYS:sim:inputSize',...
            'Input array has the wrong number of columns.')
    end
elseif isa(input,'function_handle')
    try
        testInput = input(0);
    catch
        error('POLYSYS:sim:inputSyntax', ...
            'Input function does not have appropriate calling syntax.')
    end
    if not(isequal( size(testInput), [length(sys.inputs),1] ))
        error('POLYSYS:sim:inputFunc', ...
            'Input function does not return a valid input vector.')
    end
end

nonstiffSolvers = {'ode45','ode23','ode113'};
stiffSolvers = {'ode15s','ode23s','ode23t','ode23tb'};

% Create function handles for system dynamics.
if ismember(solver,nonstiffSolvers)
    
    [stFunc,orFunc] = local_makeOdeFunction(sys,input);
    
elseif ismember(solver,stiffSolvers)
    
    [stFunc,orFunc,stJacFunc] = local_makeOdeFunction(sys,input);
    options = odeset(options,'Jacobian',stJacFunc);
    
else
    error('POLYSYS:sim:unknownSolver', ...
        'Invalid solver name ''%s''.', solver)
end

% Integrate the system dynamics.
t = cell(1,nSims);
x = cell(1,nSims);
for i = 1:nSims
    [t{i},x{i}] = feval(solver,stFunc,timeSpan,initialState(:,i),options);
end

% Compute the output (if necessary).
y = cell(1,nSims);
if nargout > 2
    for j = 1:nSims
        t_j = t{j};
        x_j = x{j};
        y_j = zeros(length(t_j),length(sys.orMap));
        for i = 1:length(t_j)
            y_j(i,:) = orFunc(t_j(i),x_j(i,:)')';
        end
        y{j} = y_j;
    end
end



%%  Convert polysys dynamics into an anonymous function.
function [f,h,Df] = local_makeOdeFunction(sys,input)

% Convert input signal to a input function.
if isa(input,'function_handle')
    inputFunc = input;
elseif size(input,1) == 1
    inputSignal = input(1,2:end)';
    inputFunc = @(tIN)( inputSignal );
else
    inputTime = input(:,1);
    inputSignal = input(:,2:end);
    inputFunc = @(tIN)( interp1(inputTime,inputSignal,tIN)' );
end

% Create a function handle for the state transition map.
stMask = pevalIndex( sys.stMap, [sys.inputs;sys.states] );
stCoef = sys.stMap.coefficient;
stDeg = sys.stMap.degmat;
nStates = size(sys.stMap,1);
f = @( tIN, xIN )( reshape( ...
        peval( stMask*[inputFunc(tIN);xIN],stCoef,stDeg),...
    nStates,1) );

% Create a function handle for the output response map.
orMask = pevalIndex( sys.orMap, [sys.inputs;sys.states] );
orCoef = sys.orMap.coefficient;
orDeg = sys.orMap.degmat;
nOutputs = size(sys.orMap,1);
h = @( tIN, xIN )( reshape( ...
        peval( orMask*[inputFunc(tIN);xIN], orCoef, orDeg ),...
    nOutputs,1) );

% If we're using a stiff solver, create a function handle for the Jacobian of f.
if nargout > 2
    jacPoly = jacobian(sys.stMap,sys.states);
    jacMask = pevalIndex(jacPoly, [sys.inputs;sys.states] );
    jacCoef = jacPoly.coefficient;
    jacDeg = jacPoly.degmat;
    Df = @(tIN,xIN)( reshape( ...
            peval( jacMask*[inputFunc(tIN);xIN],jacCoef,jacDeg),...
        nStates,nStates) );
end



%% Discrete-time local functions are commented for safe-keeping.  For
% discrete-time simulations see @polysys/dsim.m

% %% Discrete-time solver.
% function [t,x,y] = local_discSolver(sys,timeSpan,state0,input,interpMethod,nSims)
% 
% % Make sure that timeSpan is compatible with the sys.sampleTime.
% if length(timeSpan) == 1
%     
%     if sys.sampleTime == -1
%         error('POLYSYS:sim:sampleTime',...
%             'Cannot determine sample time from simulation time span.')
%     elseif sys.sampleTime > timeSpan
%         error('POLYSYS:sim:sampleTime',...
%             'Simulation time span is less than the system sampling time.')
%     else
%         tSample = sys.sampleTime;
%         timeSpan = (0:tSample:timeSpan)';
%     end
%     
% elseif length(timeSpan) == 2
%     
%     if sys.sampleTime == -1
%         if diff(timeSpan) < 0
%             error('POLYSYS:sim:sampleTime',...
%                 'Simulation time span must be increasing.')
%         else
%             tSample = diff(timeSpan);
%             timeSpan = timeSpan(:);
%         end
%     elseif sys.sampleTime > diff(timeSpan)
%         error('POLYSYS:sim:sampleTime',...
%             'Simulation time span is less than the system sampling time.')
%     else
%         tSample = sys.sampleTime;
%         timeSpan = (timeSpan(1):tSample:timeSpan(2))';
%     end
%     
% else
% 
%     % Make sure that timeSpan has evenly spaced entries.
%     diffTime = diff(timeSpan);
%     if not(all(diffTime(1)==diffTime))
%         error('POLYSYS:sim:sampleTime',...
%             'The input must be specified at regularly spaced time points.')
%     end
%     
%     if sys.sampleTime == -1 || (sys.sampleTime == diffTime(1))
%         tSample = diffTime(1);
%         timeSpan = reshape(timeSpan,length(timeSpan),1);
%     else
%         error('POLYSYS:sim:sampleTime',...
%             'System and input signal must have the same sampling time.')
%     end
% 
% end
% nSteps = length(timeSpan);
% 
% % Verify that the system input is compatible with timeSpan.
% nInputs = length(sys.inputs);
% if isempty(input)
%     u = zeros(nInputs,nSteps);
% elseif isa(input,'function_handle')
%     try
%         testInput = input(0);
%     catch
%         error('POLYSYS:sim:inputFunc', ...
%             'Input function does not have appropriate calling syntax.')
%     end
%     if isequal( size(testInput), [length(sys.input),1] )
%         u = zeros(nInputs,nSteps);
%         for k = 1:nSteps
%             u(:,k) = input(timeSpan(k));
%         end
%     else
%         error('POLYSYS:sim:inputFunc', ...
%             'Input function does not return a valid input vector.')
%     end
% elseif isnumeric(input)
%     if size(input,2) == (nInputs+1)
%         
%         if size(input,1) == 1
%             u = repmat(input(1,2:end)',1,length(timeSpan));
%         else
%             
%             if strcmp(interpMethod,'zoh')
%                 method = 'nearest';
%             else
%                 method = 'linear';
%             end
%             if (min(input(:,1)) > min(timeSpan)) || (max(input(:,1)) < max(timeSpan))
%                 error('POLYSYS:sim:input',...
%                     'Input signal not defined over the whole time span.')
%             end
%             u = interp1(input(:,1),input(:,2:end),timeSpan,method)';
%         end
%     else
%         error('POLYSYS:sim:input',...
%             'Input array has the wrong number of columns.')
%     end
% else
%     error('POLYSYS:sim:input',...
%         'System input may be empty, a double array, or a function handle.')
% end   
% 
% 
% % Convert system dynamics into functions.
% [f,g] = local_makeDiscFunction(sys);
% 
% % Run the simulation.
% t = cell(1,nSims);
% x = cell(1,nSims);
% y = cell(1,nSims);
% for i = 1:nSims
%     % Preallocate.
%     x_i = zeros(size(state0,1),nSteps);
%     y_i = zeros(size(sys.orMap,1),nSteps);
%     
%     % Load initial conditions.
%     x_i(:,1) = state0(:,i);
%     y_i(:,1) = g(x_i(:,1),u(:,1));
%     
%     % Simulate remaining time steps.
%     for j = 2:nSteps
%         x_i(:,j) = f(x_i(:,j-1),u(:,j-1));
%         y_i(:,j) = g(x_i(:,j),u(:,j));
%     end
%     
%     t{i} = timeSpan;
%     x{i} = x_i';
%     y{i} = y_i';
%     
% end
% 
% 
% function [f,g] = local_makeDiscFunction(sys)
% 
% % Create a function handle for the state transition map.
% [stInputMask,stStateMask] = pevalIndex( sys.stMap, sys.inputs, sys.states );
% stCoef = sys.stMap.coefficient;
% stDeg = sys.stMap.degmat;
% nStates = size(sys.stMap,1);
% f = @( xIN, uIN )( reshape( ...
%         peval4([stInputMask*uIN;stStateMask*xIN],stCoef,stDeg),...
%     nStates,1) );
% 
% % Create a function handle for the output response map.
% [orInputMask,orStateMask] = pevalIndex( sys.orMap, sys.inputs, sys.states );
% orCoef = sys.orMap.coefficient;
% orDeg = sys.orMap.degmat;
% nOutputs = size(sys.orMap,1);
% g = @( xIN, uIN )( reshape( ...
%         peval4([orInputMask*uIN;orStateMask*xIN],orCoef,orDeg),...
%     nOutputs,1) );
