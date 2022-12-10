function linSys = ss(sys)
%SS  Convert a linear POLYSYS model to an SS object.
%
%   LINSYS = SS(SYS) converts the POLYSYS model SYS to an SS object LINSYS.
%   The input SYS must be linear.
%
%   See also LINEARIZE, ISLINEAR, SS.

% 7.22.2007: TJW - Initial coding


if islinear(sys)

    [nOutputs,nInputs] = size(sys);
    nStates = length(sys.states);
    
    A = zeros(nStates);
    B = zeros(nStates,nInputs);
    C = zeros(nOutputs,nStates);
    D = zeros(nOutputs,nInputs);
    
    
    stCoeffs = full(sys.stMap.coefficient);
    stVarNames = sys.stMap.varname;
    orCoeffs = full(sys.orMap.coefficient);
    orVarNames = sys.orMap.varname;

    stInputCount = 0;
    orInputCount = 0;
    for i = 1:nInputs
        if ismember(sys.inputs(i).varname{1},stVarNames)
            stInputCount = stInputCount + 1;
            B(:,i) = stCoeffs(stInputCount,:)';
        end
        if ismember(sys.inputs(i).varname{1},orVarNames)
            orInputCount = orInputCount + 1;
            D(:,i) = orCoeffs(orInputCount,:)';
        end
    end
    stStateCount = 0;
    orStateCount = 0;
    for i = 1:nStates
        if ismember(sys.states(i).varname{1},stVarNames)
            stStateCount = stStateCount + 1;
            A(:,i) = stCoeffs(stInputCount+stStateCount,:)';
        end
        if ismember(sys.states(i).varname{1},orVarNames)
            orStateCount = orStateCount + 1;
            C(:,i) = orCoeffs(orInputCount+orStateCount,:)';
        end
    end
    
    if isempty(sys.sampleTime)
        linSys = ss(A,B,C,D);
    else
        linSys = ss(A,B,C,D,sys.sampleTime);
    end
    linSys.Name = sys.name;
            
else
    error('POLYSYS:nonlinear',...
        'Only linear systems can be converted to SS objects.')
end