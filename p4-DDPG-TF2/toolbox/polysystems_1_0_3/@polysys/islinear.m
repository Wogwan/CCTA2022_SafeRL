function reply = islinear(sys)
%ISLINEAR  Returns true if the polysys is linear.
%
%   REPLY = islinear(SYS) returns true (logical(1)) if the polysys object
%   SYS represents a linear system.  Otherwise, REPLY is false.

% Developer Notes:
% 7.22.2007: TJW - Initial coding.
% 7.29.2007: TJW - Fixed reply for static maps.
% 7.29.2007: TJW - Help documentation, comments.

error(nargchk(1,1,nargin,'struct'));

% Check the state transition map.
if sys.isDynamic
    stDegrees = sys.stMap.degmat;
    [stRows,stCols] = size(stDegrees);

    isSTLinear = (stRows == stCols)  &&  isequal(stDegrees,eye(stRows));
else
    isSTLinear = true(1,1);
end

% Check the output response map.
orDegrees = sys.orMap.degmat;
[orRows,orCols] = size(orDegrees);

isORLinear = (orRows == orCols)  &&  isequal(orDegrees,eye(orRows));

% Reply is true if both maps are linear.
reply = isSTLinear && isORLinear;
    