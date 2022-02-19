function sys = sys_launch(u)

pvar x1 x2;
a_rl = u;
% pvar x1 x2;
x = [x1;x2];
%% System information
sys.dom = 6;
sys.domain = [-sys.dom sys.dom -sys.dom sys.dom];
sys.dt = 0.05;
sys.g = 10;
sys.gg = [0;1];
sys.l = 1.; % Real value = 1.0
sys.m = 1.; % Real value = 1.0
sys.h_degree = 2;
sys.L_au = 6;
% Taylor based equation: taylor(10*sin(x1))
sys.f = [x2
    x1^5/12-(5*x1^3)/3+10*x1+a_rl];
% Chebyshev based equation
% sys.f = [x2
%     0.070421222500154631585012054983963*x1^5+0.00000000000000071054273576010022531416952862872*x1^4-1.6415384826007917012935521938743*x1^3-0.0000000000000021316282072803006759425085858861*x1^2+9.9876513003548918895324604250163*x1+0.0000000000000015987211554602254576530748631014+a_rl];
C1 = x1 + 1; C2 = -x1 + 1;
C3 = x2 + 8; C4 = -x2 + 8;
sys.us_region = [C1;C2;C3;C4];
sys.h = 1 - (x1^2+x1*x2+x2^2);
sys.input_limit = 8; u_out = 1;
end