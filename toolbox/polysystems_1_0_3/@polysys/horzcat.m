function sys = horzcat(sys1,varargin)
%HORZCAT  Horizontally concatenate two POLYSYS objects.
%
%   SYS = [SYS1,SYS2]  appends the inputs of SYS2 to the inputs of SYS1
%   and adds the outputs of SYS1 to the outputs of SYS2.  Thus, both 
%   systems must have the same number of outputs. This syntax can also be used 
%   to connect more than two systems (i.e. SYS = [SYS1;SYS2;SYS3;...]).
%
%   SYS = horzcat(SYS1,SYS2)  explicitly calls this method, rather than
%   using the overloaded operator syntax above.
%  
%   Note: outputs are connected according to their order in each system, so
%   that SYS.orMap(i) = SYS1.orMap(i) + SYS2.orMap(i) for i=1:"NumOutputs".
%
%   See also  polysys/vertcat, polysys/append, polysys/parallel

% 07.21.2007: TJW - Initial coding.
% 07.22.2007: TJW - Changed class name from twpds to polysys.
% 07.22.2007: TJW - Help documentation.
% 08.22.2007: TJW - Use combineNames(), matchSamplingTimes() utilities.
%                   Extended to allow more than two systems.
%                   Added nargchk, nargoutchk.


% Check the number of input and output arguments.
maxNargin = max(2,nargin);
error(nargchk(2,maxNargin,nargin,'struct'))
error(nargoutchk(0,1,nargout,'struct'))

% Make sure all items are polysys & have correct number of ouputs.
nOutputs = length(sys1.orMap);
sysArray = cell(size(varargin));
for i = 1:length(varargin)
    if isa(varargin{i},'polysys')
        sysArray{i} = varargin{i};
    else
        try
            sysArray{i} = toPolysys(varargin{i});
        catch
            error('POLYSYS:horzcat:invalidType', ...
                'Cannot convert argument %0.0f to a polysys.',i+1);
        end
    end
    if length(sysArray{i}.orMap)~=nOutputs
        error('POLYSYS:horzcat:numOutputs', ...
            'Systems must have the same number of outputs.')
    end
end

% Combine systems.
stMap = sys1.stMap;
orMap = sys1.orMap;
states = sys1.states;
inputs = sys1.inputs;
name = sys1.name;
for i = 1:length(sysArray)
    sysArray{i} = normalizeVars(sysArray{i},length(states),length(inputs));
    stMap = [stMap;sysArray{i}.stMap];
    orMap = orMap + sysArray{i}.orMap;
    states = [states;sysArray{i}.states];
    inputs = [inputs;sysArray{i}.inputs];
    if not(isempty(sysArray{i}.name))
        name = [name,sysArray{i}.name];
    end
end
sampleTime = matchSamplingTimes(sys1,sysArray{:});

sys = polysys(stMap,orMap,states,inputs,sampleTime);
sys.name = name;
