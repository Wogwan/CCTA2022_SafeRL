function varargout = size(sys,varargin)
%SIZE  Determine the size of a POLYSYS model.
%
%   D = SIZE(SYS) returns the dimensions of the POLYSYS model SYS.
%   Here, D = [NY NU] where NY is the number of outputs and NU is the
%   number of inputs.
%
%   Alternative syntax:
%   * [NY,NU] = SIZE(SYS)
%   * NY = SIZE(SYS,1)
%   * NU = SIZE(SYS,2)

% 7.21.2007: TJW - Initial coding.


%% Get system dimensions.
nOutputs = length(sys.orMap);
nInputs = length(sys.inputs);


%% Handle each input-output combination.
switch nargout
    
    case {0,1}
        
        switch nargin
            case 1
                % User wants both sizes in an array.
                varargout{1} = [nOutputs,nInputs];
            case 2
                if varargin{1} == 1
                    % User wants 'first' dimension.
                    varargout{1} = nOutputs;
                elseif varargin{1} == 2
                    % User wants 'second' dimension.
                    varargout{1} = nInputs;
                else
                    error('POLYSYS:outOfRange',...
                        'Dimension must be either 1 or 2.')
                end
            otherwise
                error('POLYSYS:nargin','Too many input arguments.')
        end
                    
    case 2
        
        % User wants both sizes.
        if nargin == 1
            varargout{1} = nOutputs;
            varargout{2} = nInputs;
        else
            error('POLYSYS:nargin',...
                'No extra inputs are needed with two output arguments.')
        end

    otherwise
        error('POLYSYS:nargout','Too many output arguments.')
end