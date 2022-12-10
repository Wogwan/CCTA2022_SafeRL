function sys = uminus(sys)
%UMINUS  Negate the outputs of a POLYSYS model.
%
%   NEWSYS = UMINUS(SYS) performs NEWSYS = -SYS, which negates the outputs
%   of the POLYSYS model SYS.
%
%   See also POLYSYS, UMINUS, TIMES.

% 7.21.2007: TJW - Initial coding.

error(nargchk(1,1,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

sys.orMap = -sys.orMap;