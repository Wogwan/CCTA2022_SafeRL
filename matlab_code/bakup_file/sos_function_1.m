function [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1(f,k,L_au,solh,V,gamma,gg)
pvar x1 x2 htol epsi;
x = [x1;x2];
%% Create corresponding decision variable
[L1,L1_Q] = sosdecvar('L1_w',monomials(x,0:L_au/2));
[L2,L2_Q] = sosdecvar('L2_w',monomials(x,0:L_au/2));
%%
[u1,u1_Q] = polydecvar('u1_w',monomials(x,0:k)); 
[u2,u2_Q] = polydecvar('u2_w',monomials(x,0:k)); 
%% CLBF
hdot = jacobian(solh, x1)*(f(1)+gg(1)*u1)+jacobian(solh, x2)*(f(2)+gg(2)*u2);
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+jacobian(V, x2)*(f(2)+gg(2)*u2);
%% Constrain :
sosconstr_1 = L1 >= 0;
sosconstr_2 = L2 >= 0;
sosconstr_3 = -Vdot>= L1*solh;
sosconstr_4 = hdot+gamma*solh-L2*solh-htol >= 0;
sosconstr_5 = htol>=0;
sosconstr = [sosconstr_1;sosconstr_2;sosconstr_3;sosconstr_4;sosconstr_5];
%% Set objection
obj = -htol;
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(sosconstr,x,obj,opts);
%% Create output
if info.feas
    SOL1 = subs(L1,dopt);
    SOL2 = subs(L2,dopt);
    SOLu1 = subs(u1,dopt);
    SOLu2 = subs(u2,dopt);
    kk = 1;
else
    SOL1 = 0;
    SOL2 = 0;
    SOLu1 = 0;
    SOLu2 = 0;
    kk = 0;
    fprintf('L1 and L2 can not find.====== ');
end
end