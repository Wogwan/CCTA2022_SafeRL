% min a_t
% s.t. h(s_{t+1})-h(s_{t}) >= 0
%      h >= 0



%% Build system
pvar x1 x2
% Taylor based equation
% f = [x2
%     -x1^5/12+(5*x1^3)/3-10*x1];
% Chebyshev based equation
f = [x2
    -0.070421222500154631585012054983963*x1^5-0.00000000000000071054273576010022531416952862872*x1^4+1.6415384826007917012935521938743*x1^3+0.0000000000000021316282072803006759425085858861*x1^2-9.9876513003548918895324604250163*x1-0.0000000000000015987211554602254576530748631014];
x = [x1;x2];
C1 = x1 + 1;
C2 = -x1 + 1;
C3 = x2 + 8;
C4 = -x2 + 8;
sys.us_region = [C1;C2;C3;C4];
%% Polynomial
[u2,uc2] = polydecvar('u_w2',monomials(x,0:2)); 
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:sys.L_au/2));
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:sys.L_au/2));
% [h1,h1_Q] = polydecvar('h1_w',monomials(x,0:sys.h_degree)); 
h = 0.01 - (x1^2+x1*x2+x2^2);
V = x1^2+x1*x2+x2^2;
Vdot = jacobian(V, x1)*(f(1))+jacobian(V, x2)*(f(2)+sys.gg(2)*u2);
hdot = jacobian(h, x1)*(f(1))+jacobian(h, x2)*(f(2)+sys.gg(2)*u2);

%% Constraint:
% pconstr_11 = L1 >= 0;
pconstr_12 = L2 >= 0;
% pconstr_21 = -Vdot - L1*h >= 0;
pconstr_22 = hdot - L2*h >= 0;
% pconstr_31 = -h
% pconstr_32 = 
% pconstr_33 = 
% pconstr_34 = 
% pconstr = [pconstr_11; pconstr_12; pconstr_21; pconstr_22];
pconstr = [pconstr_12; pconstr_22];

%% Objective
obj = sum(uc2);

%% Solver parameters
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
% opts.solver = 'sedumi';
% [info,dopt] = sosopt(pconstr,[x1;x2],opts);
[info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);

%% Check output
if info.feas
    u_output = subs(u2,dopt)
else
    fprintf('Minimum input can not find.======\n');
    return;
end