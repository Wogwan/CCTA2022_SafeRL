function varargout = polyVarsIndex(poly,varargin)

% Make sure we have two or more input arguments.
maxArgs = max(2,nargin);
error(nargchk(2,maxArgs,nargin,'struct'))

% Make sure we have fewer outputs than inputs.
error(nargoutchk(0,nargin-1,nargout,'struct'))

if isempty(poly)
    for i = 1:nargout
        varargout{i} = [];
    end
else
    polyVars = get(poly,'varname');
    for i = 1:nargout
        vars = get(varargin{i},'varname');
        vars = vars(varargin{i}.coef*(1:length(vars))');
        [commonVars,locIndex] = intersect(vars,polyVars);
        varargout{i} = locIndex;
    end
end