function sys = mtimes(sys2,sys1)
%MTIMES  Multiplication of POLYSYS models.
% 
%   SYS = MTIMES(SYS1,SYS2) performs SYS = SYS1 * SYS2. Multiplying two
%   POLYSYS models is equivalent to connecting them in series as shown below:
%
%         u ----> SYS2 ----> SYS1 ----> y 
%
%   See also POLYSYS, SERIES, TIMES.

% 7.20.2007: TJW - Initial coding.
% 7.22.2007: TJW - Changed class name from twpds to polysys.

%% Multiply two polysys objects.
if isa(sys2,'polysys')  &&  isa(sys1,'polysys')
    
    % Make sure systems have compatible dimensions.
    nOut1 = length(sys1.orMap);
    nIn2 = length(sys2.inputs);
    if nIn2 ~= nOut1
        error('POLYSYS:mtimes:innerDims',...
            'The number of inputs to the second system must match the number of outputs of the first.')
    end


    % Make sys2 states come in sequence after sys1 states.
    % (If sys1 has states x1,x2,...,xn ans sys2 has states x1,...,xm,
    % then we change sys2 states to xn+1,...,xn+m.)

    nStates1 = length(sys1.states);
    sys2 = normalizeVars(sys2,nStates1,[]);


    % Combine data together so that outputs=sys2(sys1(inputs)).

    % Combine the system mappings.
    if isempty(sys2.stMap)
        stMap = sys1.stMap;
    else
        stMap = [ sys1.stMap; robustSubs(sys2.stMap,sys2.inputs,sys1.orMap) ];
    end
    orMap = robustSubs(sys2.orMap,sys2.inputs,sys1.orMap);

    % Combine other data together.
    states = [sys1.states;sys2.states];
    inputs = sys1.inputs;
    sampleTime = matchSamplingTimes(sys1,sys2);


    % Create the new cascaded system object.
    sys = polysys(stMap,orMap,states,inputs,sampleTime);

%% Multiply polysys object and scalar.
else
    % Note: code is duplicated here, because multiplication of systems is
    % actually composition of two operators and is not commutative.
    if isa(sys1,'polysys')
        if isnumeric(sys2) && isscalar(sys2)
            sys = sys2.*sys1;
            return;
        else
            sys2 = toPolysys(sys2);
        end
    else
        if isnumeric(sys1) && isscalar(sys1)
            sys = sys1.*sys2;
            return;
        else
            sys1 = toPolysys(sys1);
        end
    end
    sys = mtimes(sys2,sys1);
end
        
    