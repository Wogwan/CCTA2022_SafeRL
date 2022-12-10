function sys = feedback(sys1,sys2,varargin)
%FEEDBACK  Feedback connection of two POLYSYS models. 
%  
%   SYS = feedback(SYS1,SYS2) computes an POLYSYS model SYS for
%   the closed-loop feedback system
%
%          u --->O---->[ SYS1 ]----+---> y
%                |                 |           y = SYS * u
%                +-----[ SYS2 ]<---+
%
%   Negative feedback is assumed and the resulting system SYS 
%   maps u to y.  To apply positive feedback, use the syntax
%   SYS = feedback(SYS1,SYS2,+1).
%
%   SYS = feedback(SYS1,SYS2,FEEDIN,FEEDOUT,SIGN) builds the more
%   general feedback interconnection:
%                      +--------+
%          v --------->|        |--------> z
%                      |  SYS1  |
%          u --->O---->|        |----+---> y
%                |     +--------+    |
%                |                   |
%                +-----[  SYS2  ]<---+
%
%   The vector FEEDIN contains indices into the input vector of SYS1
%   and specifies which inputs u are involved in the feedback loop.
%   Similarly, FEEDOUT specifies which outputs y of SYS1 are used for
%   feedback.  If SIGN=1 then positive feedback is used.  If SIGN=-1 
%   or SIGN is omitted, then negative feedback is used.  In all cases,
%   the resulting POLYSYS model SYS has the same inputs and outputs as 
%   SYS1 (with their order preserved).
%
%   See also POLYSYS/PARALLEL, POLYSYS/SERIES, POLYSYS/APPEND.

% Developer Notes:
% 08.23.2007: TJW - Initial coding & help doc.

% Check the number of input and output arguments.
error(nargchk(2,5,nargin,'struct'))
error(nargoutchk(0,1,nargout,'struct'))

% Make sure systems are both polysys.
if isa(sys1,'polysys')
    if not(isa(sys2,'polysys'))
        sys2 = toPolysys(sys2);
    end
else
    sys1 = toPolysys(sys1);
end

% Make sure systems are of the same type before we do anything else.
samplingTime = matchSamplingTimes(sys1,sys2);

% Handle each calling syntax separately.
[nOut1,nIn1] = size(sys1);
[nOut2,nIn2] = size(sys2);
switch nargin
    case 2
        feedIn = (1:nIn1)';
        feedOut = (1:nOut1)';
        sign = -1;
    case 3
        feedIn = (1:nIn1)';
        feedOut = (1:nOut1)';
        sign = varargin{1};
    case 4
        feedIn = varargin{1};
        feedOut = varargin{2};
        sign = -1;
    case 5
        feedIn = varargin{1};
        feedOut = varargin{2};
        sign = varargin{3};
end

% Make sure input arguments are valid.
[feedOut,feedIn] = checkSubscripts(sys1,feedOut,feedIn);
if length(feedIn) ~= nOut2
    error('POLYSYS:feedback:feedInDims', ...
        'Dimension mismatch on the input indices.')
end
if length(feedOut) ~= nIn2
    error('POLYSYS:feedback:feedOutDims', ...
        'Dimension mismatch on the output indices.')
end
if not( isnumeric(sign) && isscalar(sign) )
    error('POLYSYS:feedback:invalidSign', ...
        'Sign must be a numeric scalar.')
end


% Make necessary substitutions.
sys2 = normalizeVars(sys2,length(sys1.states),[]);

if not( sys1.hasDirectFeedthrough ) && not( sys2.hasDirectFeedthrough )

  w = sys1.inputs(feedIn,1) + sign*sys2.orMap;
  
  f2 = robustSubs( sys2.stMap, sys2.inputs, sys1.orMap(feedOut,1));
  f1 = robustSubs( sys1.stMap, sys1.inputs(feedIn,1), w );
  g = sys1.orMap;
  
elseif not( sys2.hasDirectFeedthrough )
  
  w = sys1.inputs(feedIn,1) + sign*sys2.orMap;
  
  f1 = robustSubs( sys1.stMap, sys1.inputs(feedIn,1), w );
  g  = robustSubs( sys1.orMap, sys1.inputs(feedIn,1), w );
  f2 = robustSubs( sys2.stMap, sys2.inputs, g(feedOut) );
 
elseif not( sys1.hasDirectFeedthrough )
  
  g2 = robustSubs( sys2.orMap, sys2.inputs, sys1.orMap(feedOut,1) );
  w = sys1.inputs(feedIn,1) + sign*g2;
  
  f1 = robustSubs( sys1.stMap, sys1.inputs(feedIn,1), w );
  f2 = robustSubs( sys2.stMap, sys2.inputs, sys1.orMap(feedOut,1) );
  g  = sys1.orMap;
  
else  
  error('POLYSYS:feedback:feedthrough', ...
        'Cannot make feedback connections when systems have direct feedthrough.')
end
  
stMap = [f1;f2];
orMap = g;
states = [sys1.states;sys2.states];
inputs = sys1.inputs;
name = combineNames(sys1,sys2);

sys = polysys(stMap,orMap,states,inputs,samplingTime);
sys.name = name;



