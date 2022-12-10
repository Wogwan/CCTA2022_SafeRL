function sys = append(sys1,varargin)
%APPEND  Group polysys models by appending their inputs and outputs.
%
%   SYS = append(SYS1,SYS2,...)  concatenates the input and output vectors
%   of SYS1 and SYS2 to form SYS.  

% 07.21.2007: TJW - Initial coding.
% 07.22.2007: TJW - Changed class name from twpds to polysys.
% 08.21.2007: TJW - Changed to allow for empty names.

% Check number of input and output arguments.
maxNargin = max(2,nargin);
error(nargchk(2,maxNargin,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

% Make sure all items are polysys
if not(isa(sys1,'polysys'))
  try
    sys1 = toPolysys(sys1);
  catch
    error('POLYSYS:append:invalidClass', ...
      'Cannot convert argument 1 to a polysys.');
  end
end
sysArray = cell(size(varargin));
for i = 1:length(varargin)
    if isa(varargin{i},'polysys')
        sysArray{i} = varargin{i};
    else
        try
            sysArray{i} = toPolysys(varargin{i});
        catch
            error('POLYSYS:append:invalidClass', ...
                'Cannot convert argument %0.0f to a polysys.',i+1);
        end
    end
end

% Make sure all system types are compatible.
samplingTime = matchSamplingTimes(sys1,sysArray{:});

% Concatenate system data to create larger system.
stMap  = sys1.stMap;
orMap  = sys1.orMap;
states = sys1.states;
inputs = sys1.inputs;
for i = 1:length(sysArray)
    if isempty(sysArray{i})
        sysArray{i} = [];
    else
    
        sysArray{i} = normalizeVars(sysArray{i},length(states),length(inputs));

        stMap  = [stMap;sysArray{i}.stMap];
        orMap  = [orMap;sysArray{i}.orMap];
        states = [states;sysArray{i}.states];
        inputs = [inputs;sysArray{i}.inputs];
    end
end

sys = polysys(stMap,orMap,states,inputs,samplingTime);
