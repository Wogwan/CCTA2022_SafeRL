function [dx,f,g] = sos_model(k)
%%SOS_MODEL records the state space systems for Chebyshev approximation and GP.
% In:
%     k      char     1 x 1 Input system index
% Out:
%     dx     syms     2 x 1 Original dynamic systems
%     f      syms     1 x 1 Polynomial term in systems
%     g      syms     1 x 1 Non-polynomial term in systems
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05

    syms x1 x2 x3;
    x = [x1 x2 x3];
    switch k
       case '2d_1'
          f = [-x(1)+x(2); 
               x(1)^2*x(2)];
          g = [0;
            1-sqrt(sqrt((exp(x(1))*cos(x(1)))^2))];
          dx = f+g;
       case '2d_2'
          f = [-x(1)^3+1*x(1)*x(2); 
               -x(1)^3*x(2)];
          g = [0;
            -cos(x(1))*sin(x(2))];
          dx = f+g;
      otherwise
          f = [ -x(1)-2*x(2)+x(1)*(4*x(1)^2-2*x(1)*x(2)+x(2)^2); 
              2*x(1)+x(2)*(4*x(1)^2-2*x(1)*x(2)+x(2)^2)];
          g = [0;
            sqrt(sqrt((exp(x(:,1))*cos(x(:,1)))^2))];
          dx = f+g;
    end
end