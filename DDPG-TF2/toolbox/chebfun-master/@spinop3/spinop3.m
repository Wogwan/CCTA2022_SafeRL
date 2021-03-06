classdef spinop3 < spinoperator
%SPINOP3   Class for representing the spatial part of 3D differential operators 
%for time-dependent PDEs.
%   SPINOP3 is a class for representing the spatial part S of a time-dependent 
%   PDE of the form u_t = S(u) = Lu + N(u) in 3D, where L is a linear 
%   operator and N is a nonlinear operator. 
%
%   S = SPINOP3(PDECHAR) creates a SPINOP3 object S defined by the string 
%   PDECHAR. Strings available include 'GL' for Ginzburg-Landau equation and 
%   'GS' for Gray-Scott equations. Other PDEs are available, see HELP/SPIN3.
%
%   S = SPINOP3(DOM, TSPAN) creates a SPINOP3 object S on DOM x TSPAN. The other
%   fields of a SPINOP3 are its linear part S.LIN, its nonlienar part S.NONLIN
%   and the initial condition S.INIT. The fields can be set via S.PROP = VALUE. 
%   See Remarks 1/2 and Example 1.
%
% Remark 1: The linear part has to be of the form 
%           
%  @(u) A*lap(u) + B*biharm(u) + C*triharm(u) + D*quadharm(u) + E*quintharm(u),
%
%           for some numbers A, B, C, D and E.
%
% Remark 2: The nonlinear part has to be of the form @(u) f(u), where f is a 
%           nonlinear nondifferential operator with constant coefficients.
%
% Example 1: To construct a SPINOP3 corresponding to the GL equation on 
%            DOM = [0 100]^3 x TSPAN = [0 70] with random initial condition, one 
%            can type
%
%            dom = [0 50 0 50 0 50]; tspan = [0 100];
%            S = spinop3(dom, tspan);
%            S.lin = @(u) lap(u);
%            S.nonlin = @(u) u - (1+1.5i)*u.*(abs(u).^2);
%            S.init = chebfun3(.1*randn(32,32,32), dom);
%
% See also SPINOPERATOR, SPIN3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop3(varargin)
            
            % Zero input argument:
            if ( nargin == 0 )
                return
            
            % One input argument:
            elseif ( nargin == 1 )
                if ( isa(varargin{1}, 'char') == 1 )
                    [L, N, dom, tspan, u0] = parseInputs(varargin{1});
                    S.lin = L;
                    S.nonlin = N;
                    S.domain = dom;
                    S.tspan = tspan;
                    S.init = u0;
                else
                    error('SPINOP3:constructor', ['When constructing a ', ...
                        'SPINOP3 with one input argument, this argument ', ...
                        'should be a STRING.'])
                end
            
            % Two input arguments:
            elseif ( nargin == 2 )
                countDouble = 0;
                for j = 1:nargin
                    item =  varargin{j};
                    if ( isa(item, 'double') == 1 && countDouble == 0 )
                        S.domain = item;
                        countDouble = 1;
                    elseif ( isa(item, 'double') == 1 && countDouble == 1 )
                        S.tspan = item;
                    else
                    error('SPINOP3:constructor', ['When constructing a ', ...
                        'SPINOP3 with two input arguments, these arguments ', ...
                        'should be DOUBLE.'])
                    end
                end
                
            % More than two input arguments:    
            else
                error('SPINOP3:constructor', 'Too many input arguments.')
            end
            
        end
        
    end
    
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    function [L, N, dom, tspan, u0] = parseInputs(pdechar)
        %PARSEINPUTS   Parse inputs when using a STRING.
        %   [L, N, DOM, TSPAN, U0] = PARSEINPUTS(PDECHAR), where PDECHAR is a 
        %   string, outputs two function handles L and N which represent the 
        %   linear and the nonlinear parts of the PDE represented by PDECHAR, 
        %   the domain DOM, the time interval TSPAN and the initial condition
        %   U0.
            
        % Ginzburg-Landau equation:
        if ( strcmpi(pdechar, 'GL') == 1 )
            L = @(u) lap(u);
            N = @(u) u - (1 + 1.5i)*u.*(abs(u).^2);
            dom = [0 50 0 50 0 50];
            tspan = [0 100];
            vals = .1*randn(32, 32, 32);
            u0 = chebfun3(vals, dom, 'trig');
            
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS') == 1 )
            L = @(u,v) [2e-5*lap(u); 1e-5*lap(v)];
            F = 0.030; K = 0.057; 
            N = @(u,v) [F*(1 - u) - u.*v.^2; -(F+K)*v + u.*v.^2];
            G = 0.5;
            dom = G*[0 1 0 1 0 1];
            tspan = [0 5000];
            u01 = @(x,y,z) 1 - exp(-500*((x-G/2.05).^2 + (y-G/2.05).^2 + ...
                (z-G/2.15).^2));
            u01 = chebfun3(u01, dom, 'trig');
            u02 = @(x,y,z) exp(-500*((x-G/2).^2 + 2*(y-G/2).^2 + 2*(z-G/2).^2));
            u02 = chebfun3(u02, dom, 'trig');
            u0 = [u01; u02];
       
        % Schnakenberg equations:
        elseif ( strcmpi(pdechar, 'SCHNAK') == 1 )
            L = @(u,v) [lap(u); 10*lap(v)];
            N = @(u,v) [3*(.1 - u + u.^2.*v); 3*(.9 - u.^2.*v)];
            G = 25;
            dom = G*[0 1 0 1 0 1];
            tspan = [0 500];
            u01 = @(x,y,z) 1 - exp(-2*((x-G/2.15).^2 + (y-G/2.15).^2 ...
                + (z-G/2.15).^2));
            u01 = chebfun3(u01, dom, 'trig');
            u02 = @(x,y,z) .9/(.1+.9)^2 + exp(-2*((x-G/2).^2 + 2*(y-G/2).^2 ...
                + 2*(z-G/2).^2));
            u02 = chebfun3(u02, dom, 'trig');
            u0 = [u01; u02];
            
        % Swift-Hohenberg equation:
        elseif ( strcmpi(pdechar, 'SH') == 1 )
            L = @(u) -2*lap(u) - biharm(u);
            N = @(u) -.9*u - u.^3;
            G = 25;
            dom = G*[0 1 0 1 0 1];
            tspan = [0 800];
            vals = .1*randn(32, 32, 32);
            u0 = chebfun3(vals, dom, 'trig');

        else
            error('SPINOP3:parseInputs', 'Unrecognized PDE.')
        end

    end