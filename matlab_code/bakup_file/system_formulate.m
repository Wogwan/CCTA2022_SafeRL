function sys = system_formulate

%% System hyperparameters
sys.dom = 6; 
sys.domain = [-sys.dom sys.dom -sys.dom sys.dom];
sys.g = 10;
sys.l = 1;
sys.m = 1;
sys.gg = [0;1/(sys.m*(sys.l^2))];
%% The first dimension
syms x1 x2
sys.f2_npoly = -sys.g/sys.l*sin(x1);
%%
pvar x1 x2
sys.V = [
    x1^2+x2^2+x1*x2
    x1^2+x2^2+x1*x2+x1^4+x2^4
    1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2
    1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1
    ];
% x1 is using the standard 1 rad
C1 = x1 + 1;
C2 = -x1 + 1;
C3 = x2 + 15;
C4 = -x2 + 15;
% sys.us_region = [C1;C2];
sys.us_region = [C1;C2;C3;C4];
end