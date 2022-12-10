function output = check_state(s1,s2,sys_bar,coe,sys_d2)

syms x1 x2
%% Dynamical term
sys_fac = [x1^3,x1^2,x1,1];
f2_sys = sys_fac.*sys_d2;
f2_sys = sum(f2_sys);
%% Input term
sys_bc_fac = [x1^2,x1*x2,x1,x2^2,x2,1];
sys_vc_fac = [x1^2,x1*x2,x2^2,1];
f2_u = sys_bc_fac.*coe;
f2_u = sum(f2_u);
%%
f = [x2
    f2_sys];

if sys_bar ~= 1
    V = sys_vc_fac.*coe;
    V = sum(V);
    P = subs(V,{x1,x2},{s1,s2});
    dV = diff(V,x1)*f(1)+diff(V,x2)*(f(2)+f2_u);
    Q = subs(dV,{x1,x2},{s1,s2});
    if P>=0 && Q<=0
        output = 1;
    else
        output = 0;
    end
else
    B = sys_bc_fac.*coe;
    B = sum(B);
    P = subs(B,{x1,x2},{s1,s2});
    dB = diff(B,x1)*f(1)+diff(B,x2)*(f(2)+f2_u);
    Q = subs(dB,{x1,x2},{s1,s2});
    if P>=0 && Q>=0
        output = 1;
    else
        output = 0;
    end
end

% dV = diff(V,x1)*f(1)+diff(V,x2)*f(2);