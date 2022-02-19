function [opt_input, sys] = s2_opt_barrier(sys)
%%
pvar x1 x2 u htol epsi;
format long
x = [x1;x2];
%%
f = s2p(sys.f);
gg = sys.gg;
V = sys.V(1);
C = sys.us_region;
dom = sys.dom;
domain = sys.domain;
C0 = sys.cv;
%%
trace_Q1 = 1; trace_Q = 0; mm = 0; kk = 1; iter = 0;
%%
sol_B = C0 - V; solh = sol_B;
%%
k_u = 2; k_h = 2; gamma = 0;
L_us = 2; 
L_au = 2;
% L_au = 6;
TRACE = []; Barrier = []; Control = [];
%%
while 1
    iter = iter+1;
    record_Q = trace_Q;
    %%
    [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_b1(f,k_u,L_au,solh,V,gamma,gg);
    if kk == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2]];
    %%
    [solh,trace_Q,kk] = sos_function_b2(iter,f,k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us);
    if kk == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
    %% Optimal the set
    kk = 1; OO = 0;
    %%
    if kk == 0
        fprintf('Advanced Barrier Function can not find.======\n');
    end
end
toc
A = [];
for iter = 1:length(Barrier)
    A = [A; [Control(iter,:) Barrier(iter)]];
end
if ~isempty(A)
    sys.perm_b = A(end,3);
    sys.b_input = A(end,2);
    opt_input = A(end,2);
    sys.exist_bar = 1;
else
    sys.exist_bar = 0;
    sys.perm_b = 0;
    sys.b_input = 0;
    opt_input = sys.c_U;
    fprintf('Using sublevel set control\n');
end
end
%%
% fprintf('Permissive B(x) is \n%s \n\n',char(vpa(p2s(A(end,3)))));
% fprintf('Control Input u1(x) is \n%s \n\n',char(vpa(p2s(A(end,1)))));
% fprintf('Control Input u2(x) is \n%s \n\n',char(vpa(p2s(A(end,2)))));