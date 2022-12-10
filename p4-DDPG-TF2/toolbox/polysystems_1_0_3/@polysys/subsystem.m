function sys = subsystem(sys,outIndex,inIndex)
%SUBSYSTEM  Get a subsystem of a POLYSYS object.
%
%   SUBSYS = subsystem(SYS,OUT,IN) keeps only the outputs and inputs
%   specified in the index arrays OUT and IN, respectively.  Any
%   unreferenced inputs are set to zero, and unreferenced outputs are
%   discarded.  OUT and IN may be vectors of logicals or positive integers.
%   The outputs and inputs of SYS my be reordered by specifying OUT and IN
%   in nonsequential order.  Also, the elements of OUT may be nonunique, 
%   but each element of IN must be unique.
%
%   Note:  subsref() calls this function directly when using the syntax
%   SUBSYS = SYS(OUT,IN).
%
%   See also POLYSYS/SUBSREF.

% Developer notes:
% 08.22.2007: TJW - Initial coding & help doc.


% Check number of input & output arguments.
error(nargchk(3,3,nargin,'struct'))
error(nargoutchk(0,1,nargout,'struct'))

[outIndex,inIndex] = checkSubscripts(sys,outIndex,inIndex);

% Keep desired portion of output response map.
if isempty(outIndex)
    error('POLYSYS:subsystem:noOutputs', ...
        'Cannot have a polysys with no outputs.')
else
    sys.orMap = sys.orMap(outIndex);
end

% Reorder & eliminate inputs as necessary.
inAll = (1:length(sys.inputs))';
nullIndex = setdiff(inAll,inIndex);
if not(isempty(nullIndex))
    nullInputs = sys.inputs(nullIndex);
    if sys.isDynamic
        sys.stMap = robustSubs(sys.stMap,nullInputs,zeros(size(nullInputs)));
    end
    sys.orMap = robustSubs(sys.orMap,nullInputs,zeros(size(nullInputs)));
end
sys.inputs = sys.inputs(inIndex);
sys = normalizeVars(sys,[],0);
