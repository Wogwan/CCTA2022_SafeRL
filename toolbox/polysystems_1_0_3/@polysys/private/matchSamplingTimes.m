function T = matchSamplingTimes(varargin)
%Private utility function for the POLYSYS class.
%
%MATCHSAMPLINGTIMES  Computes the sampling time for a group of polysys objects.
%
%   T = matchSamplingTimes(sys1,sys2,...)  compares the type
%   (continuous-time or discrete-time) of each system and then computes the
%   sampling time T for the group.

% Developer Notes:
% 08.22.2007: TJW - Total rewrite of previous algorithm to allow for arbitrary
%                   number of systems.

% Check the number of input and output arguments.
maxNargin = max(1,nargin);
error(nargchk(1,maxNargin,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

if all(cellfun('isclass',varargin,'polysys'))
    
    % Concatenate arguments into one vector.
    samplingTimeGetter = @(s)(get(s,'sampleTime'));
    Tvec = cellfun(samplingTimeGetter,varargin);
    
    % If all are continous...
    if all(Tvec==0)
        T = 0;
        
    % Else, if none are continuous...
    elseif not(any(Tvec==0))
     
        Tspecific = Tvec(Tvec>0);
        switch length(Tspecific)
            case 0
                T = -1;
            case 1
                T = Tspecific;
            otherwise
                if all(Tspecific == Tspecific(1))
                    T = Tspecific(1);
                else
                    error('POLYSYS:matchSamplingTimes:discreteMismatch',...
                        'Sampling times do not match.')
                end
        end
        
    else
        error('POLYSYS:matchSamplingTimes:systemType',...
            'Systems must all be of same type (continuous or discrete).')
    end
    
else
    error('POLYSYS:matchSamplingTimes:inputClass', ...
        'Inputs must polysys objects.')
end