function f = sos_cheb(deg,sz)
%%SOS_CHEB generates the function from Chebyshev approximation in [g]
% In:
%     deg    double   1  x  1   Input approximation truncated series degree
%     sz     double   1  x  1   Input approximation in the given square region
% Out:
%      f     syms     1 x 1     Polynomial function in systems
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

% See corresponding [g] in SOS_MODEL.m.
% We can plot the approximation result with code below.

    a = -sz; b = sz;
    y = chebfun('1-sqrt(sqrt((exp(x)*cos(x))^2))',[a,b],'splitting','on'); % Modified here with different non-polynomial g
    y_3 = minimax(y,deg); c_3 = chebcoeffs(y_3);
    % Plot the chebyshev approximation result with original function
    % plot(y,'r',y_3,'g--'); xlim([a b]);ylim([-10 10])
    syms x1;
    T = chebyshevT([0:deg],x1);
    f0 = vpa(T*c_3);
    x2 = x1/sz;
    f = subs(f0,x1,x2);
end

