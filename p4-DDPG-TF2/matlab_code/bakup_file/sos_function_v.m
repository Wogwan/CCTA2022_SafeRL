function [solu1,solu2,solL,kk]=sos_function_v(f,gg,k,k_l,V,cc)

kk = 1;
pvar x1 x2 vtol;
x = [x1;x2];
[u1,uc1] = polydecvar('u_w1',monomials(x,0:k)); 
[u2,uc2] = polydecvar('u_w2',monomials(x,0:k)); 
%%
[L,L_Q] = sosdecvar('L_w',monomials(x,0:k_l/2));
%%
Vdot = jacobian(V, x1)*(f(1)+gg(1)*u1)+ jacobian(V, x2)*(f(2)+gg(2)*u2);
%% Constraint:
pconstr_1 = L >= 0;
pconstr_2 = -Vdot-L*(cc-V) >= 0;
pconstr = [pconstr_1; pconstr_2];
%% Solve feasibility problem
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek';
[info,dopt] = sosopt(pconstr,x,opts);
%% Create output
if info.feas
    solu1 = subs(u1,dopt);
    solu2 = subs(u2,dopt);
    solL = subs(L,dopt);
else
    kk = 0;
    solu1  = 0;
    solu2  = 0;
    solL = 0;
    fprintf('Lyapunov SOS Factor L can not find.======\n');
    return;
end
end