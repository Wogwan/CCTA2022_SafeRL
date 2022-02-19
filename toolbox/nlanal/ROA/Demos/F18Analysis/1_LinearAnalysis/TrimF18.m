function [xtrim, utrim, dx] = TrimF18(x0,opts,displayon,trimcase)
%
% DESCRIPTION:
%   This file computes the Trim / Operating point for the F/A-18 aircraft dynamics
%   described in the f18full.m file.  
%
% INPUTS: 
%    
%   x0 := Trim Initial Guess for State vector
%   poly : = Set to 1 , if the polynomial model in f18full needs to
%           be trimmed, otherwise set to zero
%   displayon := 1 will result in displaying the trimmed values with cost
%   trimcase :=  Flight Condition Description 
%
% OUTPUTS: 
%    
%   xtrim := Trim value for State vector
%   utrim := Trim value for Input vector
%   dx := f(x,u) evaluated at trim values
%
% Note: 
%   u0 := Trim Initial Guess for Input vector is provided here inside the
%   function. since, we will always have them free to vary. 

% Abhijit 11/17/2009



%=======================================================
% Set Input Terms: Theye should all be free to vary; 

d_STAB0   = 0*57.2958;       % Stabilator Deflection, rad
d_RUD0    = 0*57.2958;       % Rudder Deflection, rad
d_AIL0    = 0*57.2958;       % Aileron Deflection, rad        
T0        = 14500;           % Thrust, assumed constant
u0        = [d_AIL0; d_RUD0; d_STAB0; T0]; 
IU        = [];   % All inputs are free to vary

%=======================================================
% Set which model to trim 
poly = opts.poly;
 

%=======================================================
% 1. Steady Level Flight 

if strcmpi(trimcase,'level')

    IX = [1;2;4;5;6;7];  % Hold IXth State Fixed
    [xtrim,utrim,ytrim,dx] = trim('f18trim',x0,u0,[],IX,IU,[]) ;
    LOCALdisplay(displayon,xtrim,utrim,dx)

%=======================================================
% 2. Coordinated turn 
elseif strcmpi(trimcase,'coordturn')

    IX  = [1;2;7];   % Hold IXth State Fixed
    IDX = [1:8]';    % Force IXth State Derivative to be Steady
    [xtrim,utrim,ytrim,dx] = trim('f18trim',x0,u0,[],IX,IU,[],[],IDX) ;
    LOCALdisplay(displayon,xtrim,utrim,dx)

%=======================================================
% 3. UnCoordinated turn 
else 

    IX  = [1;2;7];   % Hold IXth State Fixed
    IDX = [1:8]';    % Force IXth State Derivative to be Steady
    [xtrim,utrim,ytrim,dx] = trim('f18trim',x0,u0,[],IX,IU,[],[],IDX); 
    LOCALdisplay(displayon,xtrim,utrim,dx)

end










% This function is for display purpose

function LOCALdisplay(displayon,xtrim,utrim,dx)

if displayon

    fprintf(' \n')
    fprintf('----------------------------------------\n')
    fprintf('-------------Trim Report----------------\n')
    fprintf('----------------------------------------\n')
    fprintf(' \n')
    fprintf(' \n')
    fprintf(' Velocity (ft/s)= %f \t  dx = %f\n',xtrim(1),dx(1));
    fprintf(' sideslip (deg) = %f \t  dx = %f\n',xtrim(2)*57.2958,dx(2));
    fprintf(' aoa (deg) = %f \t \t dx = %f\n',xtrim(3)*57.2958,dx(3));
    fprintf(' p (deg/s) = %f \t \t dx = %f\n',xtrim(4)*57.2958,dx(4));
    fprintf(' q (deg/s) = %f \t \t dx = %f\n',xtrim(5)*57.2958,dx(5));
    fprintf(' r (deg/s)= %f \t \t dx = %f\n',xtrim(6)*57.2958,dx(6));
    fprintf(' phi (deg) = %f \t \t dx = %f\n',xtrim(7)*57.2958,dx(7));
    fprintf(' theta (deg) = %f \t \t dx = %f\n',xtrim(8)*57.2958,dx(8));
    fprintf(' psi (deg) = %f \t \t dx = %f\n',xtrim(9)*57.2958,dx(9));
    fprintf(' \n')
    fprintf(' \n')
    fprintf(' Thrust (lbs-ft/s^2)= %f \n',utrim(4));
    fprintf(' Aileron (deg) = %f \n',utrim(1)*57.2958);
    fprintf(' Rudder (deg) = %f  \n',utrim(2)*57.2958);
    fprintf(' Stabilator (deg) = %f \n',utrim(3)*57.2958);
   
end
