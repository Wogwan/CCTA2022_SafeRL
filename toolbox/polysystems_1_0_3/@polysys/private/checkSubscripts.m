function varargout = checkSubscripts(sys,varargin)
%Private utility function for the POLYSYS class.
%
%CHECKSUBSCRIPTS  Verifies subscripts and provides detailed error messages.
%
%   [CI1,CI2,...] = checkSubscripts(sys,I1,I2,...) verifies that the
%   indices I1,I2,... are logicals or positive integers within the bounds
%   implied by the system dimensions.  The outputs CI1,CI2,... are "clean"
%   indices, where logicals are converted directly to numerical values and
%   row vector are converted to column vectors.

% 08.22.2007: TJW - Initial coding.

% Check the number of input & output arguments.
maxArgs = ndims(sys) + 1;
error(nargchk(2,maxArgs,nargin,'struct'));
error(nargoutchk(0,nargin-1,nargout,'struct'));

% Get the size of the system. Add ones for singleton dimensions.
sysDims = size(sys);
if length(sysDims) < length(varargin)
    sysDims = [sysDims,ones(1,length(varargin)-length(sysDyims))];
end

% Check the index in each dimension.
validIndices = cell(size(varargin));
for i = 1:length(varargin);
    
    index = varargin{i};
    n = sysDims(i);
    list = (1:n)';
    
    switch class(index)
        
        case 'char'

            if strcmp(index,':')
                validIndex = list;
            else
                error('POLYSYS:checkSubscripts:invalidString', ...
                    'Subscript indices must either be real positive integers or logicals.')                
            end
            
        case 'logical'
            
            if not(isvector(index)) && not(isempty(index))
                error('POLYSYS:checkSubscripts:notVector', ...
                    'Index must be specified as a vector.')
            elseif length(index) <= n
                validIndex = list(index);
            else
                error('POLYSYS:checkSubscripts:outOfRange', ...
                    'Index exceeds system dimensions.')
            end

        otherwise

            % "Positive integers" may be many classes, but all are numeric.
            if isnumeric(index)
                if not(isvector(index)) && not(isempty(index))
                    error('POLYSYS:checkSubscripts:notVector', ...
                        'Index must be specified as a vector.')
                elseif any(index<=0) || any(index~=floor(index))
                    error('POLYSYS:checkSubscripts:invalidType', ...
                       'Subscript indices must either be real positive integers or logicals.')
                elseif any(index>n)
                    error('POLYSYS:checkSubscripts:outOfRange', ...
                    'Index exceeds system dimensions.')
                else
                    
                    validIndex = list(index);
                end
            else
                error('POLYSYS:checkSubscripts:invalidType', ...
                    'Subscript indices must either be real positive integers or logicals.')
            end

    end % End switch-case block
    
    % Restrict input indices to be unique.
    if (i==2) && (length(unique(validIndex))~=length(validIndex))
        error('POLSYS:checkSubscripts:nonuniqueInputs', ...
            'Input indices must be unique.')
    end
    
    validIndices{i} = validIndex;

end

varargout = validIndices(1:nargout);
