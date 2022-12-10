%-------------------------------------------------------------------------
% Calculation of maximum size of input given the upper bounds of the 
% size of input-output gain (gamma =3) for a 2-state system as described below.
% The example has been taken from: 
% 
% E. Summers and A. Packard. L2 gain veri?cation for interconnections of 
% locally stable systems using integral quadratic constraints. Submitted to
% CDC, 2010
%-------------------------------------------------------------------------

pvar x1 x2 u
gamma = 0.3362;   

x = [x1;x2];
f = [-x1+0.1*x1^2-0.1*x1^3+u;   -x2+u];
h = x1-x2;

Vdeg = 2;
switch Vdeg
    case 2
        zs = monomials([x;u],0:1);
        zV = monomials(x,1);
    case 4
        zs = monomials([x;u],0:2);
        zV = monomials(x,1:2);
    case 6
        zs = monomials([x;u],0:3);
        zV = monomials(x,1:3);
end

% Create L2 Gain option Object
L2gainopts = L2gainoptions(f,x,u,h,'display','on','zV',zV,'z1',zs);

% Estimate the input sixe for the induced L2 Gain given by gamma 
[R,V,s1,iter] = L2toL2gainest(f,x,u,h,gamma,L2gainopts);

