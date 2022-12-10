function rho = cheb_rho1(deg,sz)
%%CHEB_RHO generates the value [rho] to construct the Bernstein ellipse
% In:
%     deg    double   1  x  1   Input approximation truncated series degree
%     sz     double   1  x  1   Input approximation in the given square region
% Out:
%     rho    double   1  x  1    Chebyshev approximation
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

%%
% %     clear;clf
    deg = 4;                                                  % Target chebyshev truncated series degree
    sz = 2;                                                   % Target chebyshev truncated series squared region
%%
    figure(1);clf
    a = -sz; b = sz;
    x = chebfun('x',[-sz sz]);
    y = chebfun('-cos(x)^2*sin(x)',[a,b]); % Modified here with different non-polynomial g
    plot(y,'r'); hold on;
    z = exp(2i*pi*x);
    for rho = 1.3:0.01:1.36
        e = sz*(rho*z+(rho*z).^(-1))/2; 
        plot(e), hold on
    end
end