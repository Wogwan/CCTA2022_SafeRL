function sys = parallel(sys1,sys2,in1,in2,out1,out2)
%PARALLEL  Parallel interconnection of two POLYSYS models.
%
%                          +------+
%            v1 ---------->|      |----------> z1
%                          | SYS1 |
%                   u1 +-->|      |---+ y1
%                      |   +------+   |
%             u ------>+              O------> y
%                      |   +------+   |
%                   u2 +-->|      |---+ y2
%                          | SYS2 |
%            v2 ---------->|      |----------> z2
%                          +------+
%
%   SYS = parallel(SYS1,SYS2,IN1,IN2,OUT1,OUT2) connects the two POLYSYS
%   models SYS1 and SYS2 in parallel such that the inputs specified by IN1 
%   and IN2 are connected and the outputs specified by OUT1 and OUT2 are 
%   summed.  The resulting POLYSYS model SYS maps [v1;u;v2] to [z1;y;z2].  
%   The vectors IN1 and IN2 contain indexes into the input vectors of SYS1 
%   and SYS2, respectively, and define the input channels u1 and u2 in the 
%   diagram.  Similarly, the vectors OUT1 and OUT2 contain indices into 
%   the outputs of these two systems.
%
%   SYS = parallel(SYS1,SYS2) forms the standard parallel interconnection 
%   of SYS1 and SYS2 and returns SYS=SYS2+SYS1.
%
%   Note: Unlike the LTI version of PARALLEL, each input signal can only 
%   be used once.  Thus, IN1 and IN2 cannot contain repeated indices.
%   However, OUT1 and OUT2 may contain repeated indices.
%
%   See also POLYSYS/APPEND, POLYSYS/SERIES, POLYSYS/FEEDBACK.

% Developer Notes:
% 08.22.2007: TJW - Initial coding (restarted from scratch).  Help doc.

% Check the number of input & output arguments given.
error(nargchk(2,6,nargin,'struct'))
error(nargoutchk(0,1,nargout,'struct'))

% Handle various input argument combinations.
if nargin == 2
    sys = sys1+sys2;
    return;
elseif nargin == 4
    out1 = [];
    out2 = [];
elseif nargin ~= 6
    error('POLYSYS:parallel:numInputs', ...
        'Indices must be specified for both systems.')
end

% Make sure systems are compatible

% Verify the number of inputs/outputs we're trying to connect.
if length(in1) ~= length(in2)
    error('POLYSYS:parallel:dimMismatch', ...
        'Must specify the same number of inputs for each system.')
end
if length(out1) ~= length(out2)
    error('POLYSYS:parallel:dimMismatch', ...
        'Must specify the same number of outputs for each system.')
end

% Make sure indices are valid.

[out1,in1] = checkSubscripts(sys1,out1,in1);
[out2,in2] = checkSubscripts(sys2,out2,in2);

% See which inputs/outputs pass through.
[nOut1,nIn1] = size(sys1);
out1Keep = setdiff( (1:nOut1)', out1 );
in1Keep = setdiff( (1:nIn1)', in1 );
[nOut2,nIn2] = size(sys2);
out2Keep = setdiff( (1:nOut2)', out2 );
in2Keep = setdiff( (1:nIn2)', in2 );

% Create extended systems.
sys1e = subsystem(sys1,[out1Keep;out1],[in1Keep;in1]);
sys2e = subsystem(sys2,[out2;out2Keep],[in2;in2Keep]);

% Create padding so extended systems have the same dimensions.
sys1pad = toPolysys(zeros( length(out2Keep), length(in2Keep) ));
sys1pad.sampleTime = sys1.sampleTime;
sys2pad = toPolysys(zeros( length(out1Keep), length(in1Keep) ));
sys2pad.sampleTime = sys2.sampleTime;

% Finally, we added the extended & padded systems together.
sys = append(sys1e,sys1pad) + append(sys2pad,sys2e);
