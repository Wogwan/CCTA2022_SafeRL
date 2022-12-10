function sys = minus(sys1,sys2)
%MINUS  Subtract two polysys objects.
%
%   SYS = minus(SYS1,SYS2) connects the inputs and outputs of SYS1 to the 
%   corresponding inputs and outputs of uminus(SYS2). Thus, SYS1 and SYS2 must
%   have the same dimensions.
%
%   See also POLYSYS/PLUS, POLYSYS/UMINUS, POLYSYS/MTIMES, POLYSYS/TIMES.

% 07.21.2007: TJW - Initial coding.
% 07.22.2007: TJW - Changed class name from twpds to polysys.
% 08.22.2007: TJW - Use combineNames() and matchSamplingTimes() utilities.
%                 - Help documentation.

error(nargchk(2,2,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

sys = plus(sys1,uminus(sys2));
