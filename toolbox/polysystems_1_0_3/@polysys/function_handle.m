function [f,h,A,B,C,D] = function_handle(sys)
%FUNCTION_HANDLE  Covert POLYSYS object to function handles representation.
%
%   [F,G] = function_handle(SYS) converts the state transition map of SYS
%   to the function_handle F and the output response map of SYS to the
%   function_handle G.  Both F and G take three inputs:
%      xdot = F(t,x,u)
%         y = G(t,x,u)
%   where t is time x is the state at time t, and u is the input at time t.
%
%   [F,G,A,B,C,D] = function_handle(SYS) same as above but also return the
%   partial derivatives:
%      A = df/dx,  B = df/du,  C = dg/dx,  D = dg/du
%   As in the previous case, A, B, C, and D take three inputs (t,x,u).

% 6.23.2008: TJW - Initial coding.

nStates = length(sys.states);
nInputs = length(sys.inputs);
nOutputs = size(sys.orMap,1);

% Create a function handle for the state transition map.
stMask = pevalIndex( sys.stMap, [sys.inputs;sys.states] );
stCoef = sys.stMap.coefficient;
stDeg = sys.stMap.degmat;
f = @( tIN, xIN, uIN )( reshape( ...
        peval( stMask*[uIN;xIN], stCoef, stDeg),...
        nStates,1) );

% Create a function handle for the output response map.
orMask = pevalIndex( sys.orMap, [sys.inputs;sys.states] );
orCoef = sys.orMap.coefficient;
orDeg = sys.orMap.degmat;
h = @( tIN, xIN, uIN )( reshape( ...
        peval( orMask*[uIN;xIN], orCoef, orDeg ),...
        nOutputs,1) );

% Create a function handle for A = df/dx.
if nargout > 2
    APoly = jacobian(sys.stMap,sys.states);
    if isempty(APoly.var)
        ADouble = double(APoly);
        A = @(tIN,xIN,uIN)( ADouble );
    else
        AMask = pevalIndex(APoly, [sys.inputs;sys.states] );
        ACoef = APoly.coefficient;
        ADeg = APoly.degmat;
        A = @(tIN,xIN,uIN)( reshape( ...
                peval( AMask*[uIN;xIN], ACoef, ADeg ),...
                nStates,nStates) );
    end
end

% Create a function handle for B = df/du.
if nargout > 3
    BPoly = jacobian(sys.stMap,sys.inputs);
    if isempty(BPoly.var)
        BDouble = double(BPoly);
        B = @(tIN,xIN,uIN)( BDouble );
    else
        BMask = pevalIndex(BPoly, [sys.inputs;sys.states] );
        BCoef = BPoly.coefficient;
        BDeg = BPoly.degmat;
        B = @(tIN,xIN,uIN)( reshape( ...
                peval( BMask*[uIN;xIN], BCoef, BDeg ),...
                nStates,nInputs) );
    end
end

% Create a function handle for C = dg/dx.
if nargout > 4
    CPoly = jacobian(sys.orMap,sys.states);
    if isempty(CPoly.var)
        CDouble = double(CPoly);
        C = @(tIN,xIN,uIN)( CDouble );
    else
        CMask = pevalIndex(CPoly, [sys.inputs;sys.states] );
        CCoef = CPoly.coefficient;
        CDeg = CPoly.degmat;
        C = @(tIN,xIN,uIN)( reshape( ...
                peval( CMask*[uIN;xIN], CCoef, CDeg ),...
                nOutputs,nStates) );
    end
end

% Create a function handle for D = dg/du.
if nargout > 5
    DPoly = jacobian(sys.orMap,sys.inputs);
    if isempty(DPoly.var)
        DDouble = double(DPoly);
        D = @(tIN,xIN,uIN)( DDouble );
    else
        DMask = pevalIndex(DPoly, [sys.inputs;sys.states] );
        DCoef = DPoly.coefficient;
        DDeg = DPoly.degmat;
        D = @(tIN,xIN,uIN)( reshape( ...
                peval( DMask*[uIN;xIN], DCoef, DDeg ),...
                nOutputs,nInputs) );
    end
end
