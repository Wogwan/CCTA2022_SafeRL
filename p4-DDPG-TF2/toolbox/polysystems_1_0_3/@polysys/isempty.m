function reply = isempty(sys)
%ISEMPTY  True for an empty POLYSYS object.
%
%   REPLY = ISEMPTY(SYS)  returns logical true if SYS has no state
%   transition map and no output response map.

% Developer Notes:
% 07.22.2007: TJW - Initial coding, help documentation.
% 08.23.2007: TJW - Added nargchk.

error(nargchk(1,1,nargin,'struct'))

reply = isempty(sys.stMap) && isempty(sys.orMap);