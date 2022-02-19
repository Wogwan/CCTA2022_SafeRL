function sys = series(sys1,sys2,out1,in2)
%SERIES  Series interconnection of the two POLYSYS models.
%
%                                  +------+
%                           v2 --->|      |
%                  +------+        | SYS2 |-----> y2
%                  |      |------->|      |
%         u1 ----->|      |y1   u2 +------+
%                  | SYS1 |
%                  |      |---> z1
%                  +------+
%
%   SYS = series(SYS1,SYS2,OUT1,IN2) connects two POLYSYS models SYS1 and 
%   SYS2 in series such that the outputs of SYS1 specified by OUT1 are 
%   connected to the inputs of SYS2 specified by IN2.  The vectors OUT1 
%   and IN2 contain indices into the outputs and inputs of SYS1 and SYS2, 
%   respectively.  The resulting POLYSYS model SYS maps u1 to y2.
%
%   SYS = series(SYS1,SYS2)  connects SYS1 and SYS2 in cascade and returns
%      SYS = mtimes(SYS2,SYS1).
%
%   Note: Unlike the LTI version of SERIES, each input signal can only 
%   be used once.  Thus, IN2 cannot contain repeated indices.  However, 
%   OUT1 may contain repeated indices.
%
%   See also POLYSYS/APPEND, POLYSYS/PARALLEL, POLYSYS/FEEDBACK.

% Developer Notes;
% 08.23.2007: TJW - Initial coding & help doc.

% Check the number of input and output arguments.
error(nargchk(2,4,nargin,'struct'))
error(nargoutchk(0,1,nargout,'struct'))

switch nargin
    case 2
        sys = mtimes(sys2,sys1);
        return;

    case 4
        out1 = checkSubscripts(sys1,out1);
        [dummyOut2,in2] = checkSubscripts(sys2,[],in2);

        % Let mtimes() do the error checking.
        sys = subsystem(sys2,':',in2)*subsystem(sys1,out1,':');
        
    otherwise
        error('POLYSYS:series:notEnoughInputs', ...
            'Must specify a set of indices for each system.')
    
end