function newSys = set(sys,property,value)
%SET  Set properties of a POLYSYS object.
%
%   NEWSYS = SET(SYS,PROPERTY,VALUE)

% 7.20.2007: TJW - Initial coding.
% 7.22.2007: TJW - Changed class name from twpds to polysys.

%% Properties for the polysys object.

setProperties = {'name','sampleTime'};


%% Standard code for any "set" method.

if nargin == 1
    % Just display what's gettable
    disp('Settable properties:')
    for i = 1:length(setProperties)
        disp(['   ',setProperties{i}])
    end
else
    % Match property to a known property
    if ischar(property)
        [isValidProperty,location] =...
            ismember(lower(property),lower(setProperties));
    else
        error('POLYSYS:set:inputClass','Property must be a string.')
    end

    % Return the value corresponding to property
    if isValidProperty
        newSys = LOCAL_specificSet(sys,setProperties{location},value);
    else
        error('POLYSYS:set:unknownProperty',...
            'The property ''%s'' is not settable.',property)
    end
end


%% Specific code for polysys object.

function sys = LOCAL_specificSet(sys,property,value)

switch property

    case 'name'
        
        % Check the system name for consistency
        if not(ischar(value))  ||  size(value,1) > 1
            error('System name must be a string.')
        end
        sys.name = value;
        
        
    case 'sampleTime'
        if isempty(value)
            value = 0;
        end
        if not(isscalar(value)) || ( (value<0) && (value~=-1) )
            error('System sample time must be a positive scalar.')
        else
            sys.sampleTime = value;
        end
        
    otherwise
        error(['The property ''',property,''' does not exist.'])

end % End switch-case