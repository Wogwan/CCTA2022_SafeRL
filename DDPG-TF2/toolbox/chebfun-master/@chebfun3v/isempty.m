function out = isempty(F)
%ISEMPTY   Empty Boolean check for a CHEBFUN3V object.
%   ISEMPTY(F) returns 1 if every component of F is an empty CHEBFUN3V, and
%   return 0 otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F.components) )
    out = 1; 
    return
end

% Take isempty of each component:
out = cellfun(@isempty, F.components, 'UniformOutput', false);
out = all(cell2mat(out));

end