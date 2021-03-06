function h = ldivide(f, g)
%.\   Pointwise left array divide for CHEBFUN3T objects.
%   F.\G returns (1/F)*G, if G is a CHEBFUN3T and F is a double.
%
%   F.\G returns G/F, if G is double and F is a CHEBFUN3T, but this does
%   not work if F becomes numerically close to zero.
% 
%   F.\G is the same as the command ldivide(F, G).
%
% See also CHEBFUN3T/RDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = rdivide(g, f);

end