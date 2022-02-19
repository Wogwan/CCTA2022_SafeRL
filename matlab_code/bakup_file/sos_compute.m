function y = sos_compute(u,input_1,input_2)
u = 0.1; input_1 = 0.1; input_2 = 0.1;
%%
coder.extrinsic('sym');
coder.extrinsic('syms');
coder.extrinsic('pvar');
coder.extrinsic('mpower');
coder.extrinsic('polynomial');
coder.extrinsic('sos_poly_add');
coder.extrinsic('sos_fixed_barrier');
coder.extrinsic('Try');
coder.extrinsic('catch');
coder.extrinsic('p2s');
coder.extrinsic('subs');
%%
u_store=polynomial; u_out=polynomial; u_sos=syms;
%% Initial
% Polynomial
a_x1=polynomial;
x1=polynomial; x2=polynomial; f1=polynomial; 
% Struct and double
sys.g=0; sys.gg=0; sys.l=0; sys.m=0;
sys.h_degree=0; sys.L_au=0;
sys.input_limit=0; u_out=0;
pvar x1 x2;
a_rl = u;
% pvar x1 x2;
%% System information
sys.g = 10; sys.gg = [0;1];
sys.l = 1.; % Real value = 1.0
sys.m = 1.; % Real value = 1.0
sys.h_degree = 2; sys.L_au = 6;
sys.input_limit = 8; u_out = 1;
% Taylor based equation: taylor(10*sin(x1))
f1 = x2;
% f21=polynomial; f22=polynomial; f23=polynomial; 
% f21 = a_x1^5/12; f22 = -5/3*a_x1^3; f23=10*a_x1;
f2 = polynomial;
f2 = sos_poly_add(a_rl);
% f2 = f21+f22+f23+a_rl;
% Chebyshev based equation
% sys.f = [x2
%     0.070421222500154631585012054983963*x1^5+0.00000000000000071054273576010022531416952862872*x1^4-1.6415384826007917012935521938743*x1^3-0.0000000000000021316282072803006759425085858861*x1^2+9.9876513003548918895324604250163*x1+0.0000000000000015987211554602254576530748631014+a_rl];
% C1 = x1 + 1; C2 = -x1 + 1;
% C3 = x2 + 8; C4 = -x2 + 8;
% sys.us_region = [C1;C2;C3;C4];
% sys.h = 1 - (x1^2+x1*x2+x2^2);

while 1
    sys.input_limit = sys.input_limit - 0.1;
    sys.input_limit
    %% Find optimal L1
    L1 = sos_feasu(sys,f2);
    %% Find minimize U
    u_out = sos_min_u_sosp(sys, L1, f2);
    u_store = u_out;
    if size(u_out) == size(u_store)
        continue
    else
        break
    end
    %     try
    %         if class(u_out) == 'polynomial'
    %             u_store = u_out;
    %             continue
    %         end
    %     catch
    %         u_store
    %         sys.input_limit
    %         break
    %     end
end
syms x1 x2
u_sos = p2s(u_store);
y = subs(u_sos, {x1,x2}, {input_1,input_2});
end