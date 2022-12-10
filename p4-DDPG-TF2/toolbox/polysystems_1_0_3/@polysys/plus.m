function sys = plus(sys1,sys2)
%PLUS  Add two polysys objects together.
%
%   SYS = plus(SYS1,SYS2) connects the inputs and outputs of SYS1 and SYS2
%   together to form a system of the same size.  Thus, SYS1 and SYS2 must
%   have the same dimensions.
%
%   See also POLYSYS/MINUS, POLYSYS/UMINUS, POLYSYS/MTIMES, POLYSYS/TIMES.

% 07.21.2007: TJW - Initial coding.
% 07.22.2007: TJW - Changed class name from twpds to polysys.
% 08.22.2007: TJW - Use combineNames() and matchSamplingTimes() utilities.
%                 - Help documentation.

% Verify the number of input and output arguments.
error(nargchk(2,2,nargin,'struct'))
error(nargoutchk(0,1,nargout,'struct'))

% Make sure both systems are polysys.
if isa(sys1,'polysys')
    if not(isa(sys2,'polysys'))
        try
            sys2 = toPolysys(sys2);
        catch
            error('POLYSYS:plus:invalidClass', ...
                'Cannot convert RHS argument to a polysys.')
        end
    end
elseif not(isa(sys1,'polysys'))
    try
        sys1 = toPolysys(sys1);
    catch
        error('POLYSYS:plus:invalidClass', ...
            'Cannot convert LHS argument to a polysys.')
    end
end

% Verify that these systems have the same dimensions.
if not(isequal( size(sys1), size(sys2) ))
    error('POLYSYS:plus:dimensions', ...
        'Systems must be of the same size.')
end

% Prepare sys2 for interconnection (avoid naming conflicts).
sys2 = normalizeVars(sys2,length(sys1.states),[]);

% Concatenate dynamics and add outputs.
stMap = [ sys1.stMap; sys2.stMap ];
orMap = sys1.orMap + sys2.orMap;
states = [ sys1.states; sys2.states ];
inputs = sys1.inputs;
name = combineNames(sys1,sys2);
sampleTime = matchSamplingTimes(sys1,sys2);
sys = polysys(stMap,orMap,states,inputs,sampleTime);
sys.name = name;
