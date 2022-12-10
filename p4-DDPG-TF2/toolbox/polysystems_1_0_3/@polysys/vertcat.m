function sys = vertcat(sys1,varargin)
%VERTCAT  Vertically concatenate two POLYSYS objects.
%
%   SYS = [SYS1;SYS2]  appends the outputs of SYS2 to the outputs of SYS1
%   and 'connects' the inputs of SYS1 and SYS2.  Thus, both systems must
%   have the same number of inputs.  This syntax can also be used to
%   connect more than two systems (i.e. SYS = [SYS1;SYS2;SYS3;...]).
%
%   SYS = vertcat(SYS1,SYS2,...)  explicitly calls this method, rather than
%   using the overloaded operator syntax above.
%  
%   Note: inputs are connected according to their order in each system, so
%   that SYS1.inputs(i) = SYS2.inputs(i) for i = 1:NumInputs.
%
%   See also  polysys/horzcat.m, polysys/append.m, polysys/parallel.m

% 07.21.2007: TJW - Initial coding.
% 07.22.2007: TJW - Changed class name from twpds to polysys.
% 07.22.2007: TJW - Help documentation.
% 08.22.2007: TJW - Use combineNames(), matchSamplingTimes() utilities.
%                   Extended to allow more than two systems.


% Verify inputs.
if nargin < 2
    error('POLYSYS:vertcat:oneInput', ...
        'Must specify at least two systems for vertcat.')
end


% Make sure all items are polysys & have correct number of inputs.
nInputs = length(sys1.inputs);
sysArray = cell(size(varargin));
for i = 1:length(varargin)
    if isa(varargin{i},'polysys')
        sysArray{i} = varargin{i};
    else
        try
            sysArray{i} = toPolysys(varargin{i});
        catch
            error('POLYSYS:vertcat:invalidType', ...
                'Cannot convert argument %0.0f to a polysys.',i+1);
        end
    end
    if length(sysArray{i}.inputs)~=nInputs
        error('POLYSYS:vertcat:numInputs', ...
            'Systems must have the same number of inputs.')
    end
end

% Combine systems.
stMap = sys1.stMap;
orMap = sys1.orMap;
states = sys1.states;
inputs = sys1.inputs;
name = sys1.name;
for i = 1:length(sysArray)
    sysArray{i} = normalizeVars(sysArray{i},length(states),[]);
    stMap = [stMap;sysArray{i}.stMap];
    orMap = [orMap;sysArray{i}.orMap];
    states = [states;sysArray{i}.states];
    if not(isempty(sysArray{i}.name))
        name = [name,sysArray{i}.name];
    end
end
sampleTime = matchSamplingTimes(sys1,sysArray{:});

sys = polysys(stMap,orMap,states,inputs,sampleTime);
sys.name = name;
