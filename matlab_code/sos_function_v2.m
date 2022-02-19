function [cc,kk,solu1,solu2]=sos_function_v2(f,gg,k,k_l,V,C,dom,solL,figure_id)
%%
kk = 1;
domain = [-dom dom -dom dom];
pvar x1 x2 cc;
x = [x1;x2];
%%
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:k_l/2)); 
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:k_l/2));
[L3,L3_Q] = sosdecvar('L3_w',monomials(x,0:k_l/2)); 
[L4,L4_Q] = sosdecvar('L4_w',monomials(x,0:k_l/2)); 
%%
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k)); 
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k));
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
%% Constraint:
pconstr_21 = L1 >= 0;
pconstr_22 = L2 >= 0;
pconstr_1 = -Vdot-solL*(cc-V) >= 0;
pconstr_31 = -(cc-V)+C(1)*L1 >= 0;
pconstr_32 = -(cc-V)+C(2)*L2 >= 0;
pconstr_4 = cc >= 0;
%%
input_bundary = 4;
pconstr_u = [];
for i = 1:length(uc2)
    pconstr_u = [pconstr_u; uc2(i)<=input_bundary; uc2(i)>=-input_bundary];
end
%%
if length(C) == 2
    pconstr = [pconstr_21;pconstr_22;pconstr_1;pconstr_31;pconstr_32;pconstr_4];
%     pconstr = [pconstr_u;pconstr_21;pconstr_22;pconstr_1;pconstr_31;pconstr_32;pconstr_4];
elseif length(C) == 3
    pconstr_33 = -(cc-V)+C(3)*L3 >= 0;
    pconstr_23 = L3 >= 0;
    pconstr = [pconstr_21;pconstr_22;pconstr_23;pconstr_1;pconstr_31;pconstr_32;pconstr_33;pconstr_4];
%     pconstr = [pconstr_u; pconstr_21;pconstr_22;pconstr_23;pconstr_1;pconstr_31;pconstr_32;pconstr_33;pconstr_4];
elseif length(C) == 4
    pconstr_33 = -(cc-V)+C(3)*L3 >= 0;
    pconstr_23 = L3 >= 0;
    pconstr_34 = -(cc-V)+C(4)*L4 >= 0;
    pconstr_24 = L4 >= 0;
    pconstr = [pconstr_21;pconstr_22;pconstr_23;pconstr_24;pconstr_1;pconstr_31;pconstr_32;pconstr_33;pconstr_34;pconstr_4];
%     pconstr = [pconstr_u; pconstr_21;pconstr_22;pconstr_23;pconstr_24;pconstr_1;pconstr_31;pconstr_32;pconstr_33;pconstr_34;pconstr_4];
else
    fprintf('Constraints vector does not match.======\n');
end
%% Set objection
obj = -cc;
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
% opts.solver = 'sedumi';
% opts.solver = 'sdpam'; opts.solveropts.epsilonStar = 1e-9;
[info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
%%
if info.feas
    cc = subs(cc,dopt);
    solu1 = subs(u1,dopt);
    solu2 = subs(u2,dopt);
%     refreshdata; drawnow;
else
    kk = 0;
    cc = 0;
    solu1 = 0;
    solu2 = 0;
%     fprintf('Suitable sublevel set can not find.======\n');
    return;
end
end