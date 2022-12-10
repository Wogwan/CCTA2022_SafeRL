function [A, x_num] = cheb_max(f,sz)
%%CHEB_MAX returns the maximum value and index from Cheb-Appro function
% In:
%     f      syms     1 x 1     Polynomial function in systems
%     sz     double   1  x  1   Input approximation in the given square region
% Out:  
%     A      double   1  x  1   Output abs(maximum) value of the given [f]
%     x_num  double   1  x  1   Ouput maximum value index in the given domain
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

    syms x1 x;
%     sz = 2;
%     f = 1 - sqrt(abs((exp(x1)*cos(x1))));
    x11 = linspace(-sz,sz,1000);
    K = double(subs(f,x1,x11))';
    A = max(K);
    B = min(K);
    if abs(A)>=abs(K)
        num = find(K==A);
        x_num = x11(num);
    else
        num = find(K==B);
        x_num = x11(num);
        A = B;
    end
end