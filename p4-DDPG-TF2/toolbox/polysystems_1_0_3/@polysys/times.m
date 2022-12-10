function sys = times(A,B)
%TIMES  Multiply a POLYSYS model by a scalar.
%
%   NEWSYS = M*SYS  multiplies the outputs of the POLYSYS model SYS by the
%   scalar M.
%
%   NEWSYS = SYS*M  multiples the inputs of the POLYSYS model SYS by the
%   scalar M.
%
%   See also POLYSYS, MTIMES.


% 5.14.2009: TJW - Initial coding.

%% Determine which input is a multiplier.

if isa(A,'polysys')
    sys = A;
    multiplier = B;
    onLeft = false(1);
else
    sys = B;
    multiplier = A;
    onLeft = true(1);
end

%% Verify input.

if not(isnumeric( multiplier ))
  error('POLYSYS:times:inputClass',...
    'Multiplier must be a numeric scalar.')
end
if not(isscalar( multiplier ))
  error('POLYSYS:times:nonscalarInput',...
    'Element-by-element multiplication requires a scalar multiplier.');
end

%% Perform multiplication.

if onLeft
  
  sys.orMap = multiplier * sys.orMap;
  
else

  % Create a static system that scales the inputs.
  numInputs = length( sys.inputs );
  multsys = multiplier * polysys( eye( numInputs ) );

  % Use series interconnection to affection inputs.
  sys = series( multsys, sys );
  
end
