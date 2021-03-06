function out = isempty(f) 
%ISEMPTY   True for empty CHEBFUN3T objects.
%   ISEMPTY(F) returns 1 if F is an empty CHEBFUN3T object and 0 otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isempty(f.coeffs) && isempty(f.domain) && isempty(f.vscale); 

end