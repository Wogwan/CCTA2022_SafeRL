function varargout = pevalIndex(poly,varargin)
%Private utility function for the POLYSYS class.
%
%PEVALINDEX  Determines which variables are being used in a POLYNOMIAL object.
%
%   T = pevalIndex(P,V) determines which variables in V are being used in 
%   polynomial P.  Here the polynomial object V is a N-by-1 vector of monomials.
%   Suppose that M variables in V are used in P.  Then, the output T is a 
%   M-by-N matrix such that T*V gives the appropriate subset of variables.
%
%   [T1,T2,...] = pevalIndex(P,V1,V2,...) is the same as above where Ti*Vi
%   gives the subset of Vi that is used in P.

% 08.22.2007: TJW - Initial coding.
% 12.18.2007: TJW - Added comments and help info.

% Get the names of the variables used in poly.
polyNames = poly.varname;


for i = 1:length(varargin)

    % Determine which elements of vars are used in poly.
    vars = varargin{i};
    varNames = vars.varname;
    logicalIndex = ismember(varNames,polyNames);

    % Convert the logical index into subscript indices.
    nVars = length(vars);
    varList = (1:nVars)';
    index = varList(logicalIndex);
    
    % Return a matrix M, so that M*vars = {vars used in poly}.
    n = length(index);
    varargout{i} = full(sparse( (1:n)',index,ones(n,1),n,nVars));
    % The following command is faster, but fails when n==0.
%     varargout{i} = accumarray( [(1:n)',index], ones(n,1), [n,nVars] );
   
end
