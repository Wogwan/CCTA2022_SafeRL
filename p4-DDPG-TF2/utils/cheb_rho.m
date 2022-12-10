function rho = cheb_rho(deg,sz,F)
%%CHEB_RHO generates the value [rho] to construct the Bernstein ellipse
% In:
%     deg    double   1  x  1   Input approximation truncated series degree
%     sz     double   1  x  1   Input approximation in the given square region
% Out:
%     rho    double   1  x  1    Chebyshev approximation
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

%%
%     clear;clf
%     deg = 4;                                                  % Target chebyshev truncated series degree
%     sz = 2;                                                      % Target chebyshev truncated series squared region
%%
    figure(2);clf
    a = -sz; b = sz;
    x = chebfun('x',[-sz sz]);
    y = chebfun(F,[a,b]); % Modified here with different non-polynomial g
    [y_deg, ~] = minimax(y,deg); 
    plot(y,'b'); hold on;
    plot(y_deg,'g');
    z = exp(2i*pi*x);
    for rho = 3.4:0.01:3.5
        e = sz*(rho*z+(rho*z).^(-1))/2; 
        plot(e), hold on
    end
end