function value = get(sys,property)
%GET  Access the properties of a POLYSYS model.
%
%   VALUE = get(SYS,PROP)  Returns the appropriate value if PROP is a 
%   publicly accessible property of SYS.  Type get(SYS) to see which
%   properties can be accessed in this manner.  
%
%   See also POLYSYS/SET

% Developer Notes:
% 07.19.2007: TJW - Initial coding.
% 07.22.2007: TJW - Changed class name from twpds to polysys.
% 08.23.2007: TJW - Help documentation.

%% Check the number of input and output arguments.

error(nargchk(1,2,nargin,'struct'))
error(nargoutchk(0,1,nargout,'struct'))


%% Properties for the polysys object.

getProperties = {'stMap','orMap','states','inputs',...
                   'name','sampleTime','isDynamic','hasDirectFeedthrough',...
                   'stStateIndex','stInputIndex','orStateIndex','orInputIndex'};

               
%% Standard code for any "get" method             

if nargin == 1
    % Just display what's gettable
    disp('Gettable properties:')
    for i = 1:length(getProperties)
        disp(['   ',getProperties{i}])
    end
else
    % Match property to a known property
    if ischar(property)
        [isValidProperty,location] =...
            ismember(lower(property),lower(getProperties));
    else
        error('POLYSYS:get:inputClass','Property must be a string.')
    end

    % Return the value corresponding to property
    if isValidProperty
        value = sys.(getProperties{location});
    else
        error('POLYSYS:get:unknownProperty',...
            'The property ''%s'' is not gettable.',property)
    end
end
