function varargout = linearize(sys,x,u)
%LINEARIZE  Linearize a POLYSYS object at a given operating point.
%
%   LINSYS = linearize(SYS,X,U) linearizes the POLYSYS object SYS at the
%   point (X,U) and returns the result as the SS object LINSYS.  Note that
%   LINSYS will have the same sampling time as SYS.
%
%   [A,B,C,D] = linearize(SYS,X,U) returns the matrix data that describes
%   the linear system.
%
%   [...] = linearize(SYS) same as above except that X and U are taken to
%   be vectors of zeros.

% 2008-06-25: TJW - Initial coding.
% 2008-06-17: TJW - Added variable input and output syntax.


% Use zeros if no operating point is specified.
if nargin < 3
  u = zeros(size(sys.inputs));
end
if nargin < 2
  x = zeros(size(sys.states));
end

% Verify that the inputs have the correct dimensions.
if not(isequal(size(x),size(sys.states)))
  error('POLYSYS:linearize:inputSize',...
    'State vector has wrong dimensions.')
end
if not(isequal(size(u),size(sys.inputs)))
  error('POLYSYS:linearize:inputSize',...
    'State vector has wrong dimensions.')
end

% Get the function handle representation of the system.
[fFun,gFun,AFun,BFun,CFun,DFun] = function_handle(sys);


% Evaluate functions at the given (x,u) point.
A = AFun(0,x,u);
B = BFun(0,x,u);
C = CFun(0,x,u);
D = DFun(0,x,u);

% Create the relevant outputs.
if nargout < 2
  % Return ss object if only one output is needed.
  varargout{1} = ss(A,B,C,D,sys.sampleTime);
else
  % Return all four matrices if more than one output is requested.
  varargout = {A,B,C,D};
end
