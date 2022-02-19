pvar x1 x2
%% System hyperparameters
sys.dom = 6; 
sys.domain = [-sys.dom sys.dom -sys.dom sys.dom];
sys.g = 10;
sys.l = 1;
sys.m = 1;
sys.gg = [0;1/(sys.m*(sys.l^2))];
%% The first dimension
sys.f2_npoly = -sys.g/sys.l*sin(x1);
%%
pvar x1 x2
sys.V = [
    x1^2+2*x2^2+x1*x2
    x1^2+x2^2+x1*x2+x1^4+x2^4
    1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2
    1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1
    ];
sys.sub_level = 4;
% Unsafe region
C1 = x1 + 1;
C2 = -x1 + 1;
C3 = x2 + 4;
C4 = -x2 + 1;
% sys.us_region = [C1;C2];
sys.us_region = [C1;C2;C3;C4];