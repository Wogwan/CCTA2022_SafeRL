function [opt_input, sys] = s3_opt_lyapunov(sys)
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
dom = sys.dom; domain = sys.domain;
f = s2p(sys.f);
gg = sys.gg;
C = sys.us_region;
V0 = sys.V(1); CC0 = sys.cv;
u1 = 0; u2 = sys.b_input;
B = sys.perm_b;
%%
c0 = 1; cc = 1.1; epsi = 1e-6; figure_id = 111;
%%
% figure(figure_id);clf;hold on;
% hold on; [~,~]=pcontour(B,0,domain,'g');
% [~,~]=pcontour(V0,double(CC0),domain,'m');
%% Hyperparameters of the SOSP @ CBF -> V
V_us = 4; V_au = 4; V_degree = 4; gamma = 0; k_u_V = 2; k_l_au = 4;
%%
[V, kk] = sos_function_opt_v1(f,gg,B,u1,u2,V_au,V_us,V_degree,C,gamma);
%%
if kk == 0
    fprintf('Suitable Lyapunov function can not find.======\n');
end
%% TEST FOR Sublevel Set
[a1,b1] = coeffs(p2s(V));
ccc = double(vpa(a1(end)));
C0 = ccc; cc = ccc+1;
solU = []; v_c = []; iter = 0;
%%
% figure(figure_id);hold on;
% [~,~]=pcontour(V,double(C0),domain,'r');
%%
while double(cc)-double(C0) >= epsi
    iter = iter + 1;
    if iter ~= 0
        C0 = cc;
    end
    [solL,kk]= sos_function_opt_v2(f,gg,k_u_V,k_l_au,V,cc);
    if kk == 0
        break
    end
    [cc,kk,solu1,solu2] = sos_function_opt_v3(f,gg,k_u_V,k_l_au,V,C,dom,solL,ccc,figure_id);
    v_c = [v_c; double(cc)];
    solU = [solU;[solu1,solu2]];
    if kk == 0
        %             figure(figure_id);hold on;
        %             [~,~]=pcontour(V,v_c(end),domain,'b');
        break
    end
end
%% Start to compute the control barrier function
sys.V_opt = V;
sys.cv_opt = v_c(end);
sys.initB_opt = sys.cv_opt - sys.V_opt;
sys.V_opt_input = solU(end,2);
opt_input = solU(end,2);
% sol_B = v_c(end) - N_Lya(end);
% hold on;[~,~]=pcontour(sol_B,0,domain,'k');
end
%%
% fprintf('Sublevel set of V∗(x) is \n%s \n\n',char(vpa(p2s(V))));
% fprintf('Sublevel set of V∗(x) is %.14f. \n',v_c(end));