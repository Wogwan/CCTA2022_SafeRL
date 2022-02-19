function sys = blkdiag(varargin)
%BLKDIAG  Group polysys models by appending their inputs and outputs.
%
%   This function is the same as append.
%
%   See also APPEND.

% Developer Notes:
% 08.21.2007: TJW - Initial coding.

sys = append(varargin{:});
