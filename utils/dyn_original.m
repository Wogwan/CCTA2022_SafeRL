function dxdt = dyn_original(~,x)
%DYN2D time invariant 2D nonlinear time continous dynamics
% In:
%    t      ~       time (not used)
%    x    2 x N     state
% Out:
%    dxdt 2 x N     state derivative
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

% From main_total.m 25 line: dotx1 chebyshev approximation result with 4th order remainder
% dxdt(1,:) = x(2,:)-1.0*x(1,:);
% dxdt(2,:) = x(1,:)^2*x(2,:) - 0.38410921204577513909939057157317*x(1,:) + 1.2665468516202070398124490679947*x(1,:)^2 + 0.058711157666221952533547323582752*x(1,:)^3 - 0.28781069452828345056616399233462*x(1,:)^4 - 0.29722062523126249744542803910008;

% From main_total.m 25 line: dotx1 chebyshev approximation result with 6th order remainder
% dxdt(1,:) = x(2,:)-1.0*x(1,:);
% dxdt(2,:) = x(1,:)^2*x(2,:)-0.38410921204577513909939057157317*x(1,:)+1.2665468516202070398124490679947*x(1,:)^2+0.058711157666221952533547323582752*x(1,:)^3-0.28781069452828345056616399233462*x(1,:)^4-0.14416490303838649933432236593944;

dxdt(1,:) = x(2,:)-1.0*x(1,:);
dxdt(2,:) = x(1,:)^2*x(2,:)+1-sqrt(sqrt((exp(x(1,:))*cos(x(1,:)))^2));

end