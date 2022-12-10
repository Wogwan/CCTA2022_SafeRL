function sys = normalizeVars(sys,sOffset,iOffset)
%NORMALIZEVARS  Rename variables of a polysys object in a standardized way.
%
%   NORMSYS = normalizeVars(SYS,N,M)  renames the states in SYS to
%   (xN,xN+1,...) and the inputs to (uM,uM+1,...).  The ordering of each
%   vector is preserved.  If N or M is empty, the variables will be left
%   unchanged.

% Developer Notes:
% 08.21.2007: TJW - Combined setInputIndices & setStateIndices. Changed
%                   code to create new variable vectors directly with
%                   polynomial().

% Verify the number of input/output arguments.
error(nargchk(2,3,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

if nargin < 3
    iOffset = [];
end

normalizeStates = not(isempty(sOffset)) && not(isempty(sys.states));
normalizeInputs = not(isempty(iOffset)) && not(isempty(sys.inputs));

% Set state indices.
if normalizeStates

    nStates = length(sys.states);
    
    % Create vector with new state names.
    newStateNames = cell(nStates,1);
    for i = 1:nStates
        newStateNames{i,1} = sprintf('x%0.0f',i+sOffset);
    end
    newStates = polynomial(eye(nStates),eye(nStates),newStateNames,[nStates,1]);
    
else
    newStates = [];
end

% Set input indices.
if normalizeInputs
    
    nIn = length(sys.inputs);
   
    % Create vector with new state names.
    newInputNames = cell(nIn,1);
    for i = 1:nIn
        newInputNames{i,1} = sprintf('u%0.0f',i+iOffset);
    end
    newInputs = polynomial(eye(nIn),eye(nIn),newInputNames,[nIn,1]);

else
    newInputs = [];
end

% Do substitutions.
if normalizeStates
    
%     % Here's an idea of how the state indices might streamline this routine:
%     stVarNames = sys.stMap.varname;
%     stVarNames(end-length(sys.stStateIndex)+1:end,1) = ...
%         newStateNames(sys.stStateIndex);
%     sys.stMap.varname = stVarNames;
    
    
    % Make substitutions where necessary.
    if sys.isDynamic
        sys.stMap = robustSubs(sys.stMap,sys.states,newStates);
    end
    sys.orMap = robustSubs(sys.orMap,sys.states,newStates);
    sys.states = newStates;
end
if normalizeInputs
    % Make substitutions where necessary.
    if sys.isDynamic
        sys.stMap = robustSubs(sys.stMap,sys.inputs,newInputs);
    end
    sys.orMap = robustSubs(sys.orMap,sys.inputs,newInputs);
    sys.inputs = newInputs;
end
