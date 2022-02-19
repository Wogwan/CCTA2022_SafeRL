function L1 = sos_feasu(sys,f2)
coder.extrinsic('sym');
coder.extrinsic('syms');
coder.extrinsic('pvar');
coder.extrinsic('struct');
coder.extrinsic('nvar');
coder.extrinsic('monomials');
coder.extrinsic('polynomial');
coder.extrinsic('polydecvar');
coder.extrinsic('sosdecvar');
coder.extrinsic('subsref');
coder.extrinsic('jacobian');
coder.extrinsic('sosopt');
coder.extrinsic('subs');
coder.extrinsic('sos_compute');
coder.extrinsic('sos_poly_add');
coder.extrinsic('sos_fixed_barrier');
coder.extrinsic('sos_calculate_hdot2');
coder.extrinsic('sos_constrain_2');
%%
pvar x1 x2
x2=polynomial; 
f1 = x2;
%%
f1=polynomial;
x1=polynomial; u2=polynomial; uc2=polynomial; 
hdot1=polynomial; hdot2=polynomial; 
pconstr=polynomial; pconstr_12=polynomial; pconstr_22=polynomial; 
info.feas=0; dopt=struct; opts.from=''; opts.solve='';
u_degree=2; L_au=6;
h=polynomial;
% L1=0; L1_Q=0; 
% f21 = f2;
% f22 = sos_poly_add(a_rl);
% a_x1=pvar('x1'); x1=a_x1
% a_x2=pvar('x2'); x2=a_x2;
%% Polynomial
u2 = polynomial;
[u2,uc2] = polydecvar('u_w2',monomials([x1;x2],0:u_degree)); 
L1 = polynomial;
[L1,L1_Q] = sosdecvar('L1_w',monomials([x1;x2],0:L_au/2));
h = sos_fixed_barrier;
hdot = x2;
hdot = sos_calculate_hdot2(h,f1,f2,sys.gg,u2);
% hdot2 = jacobian(h, x2)*(f21+sys.gg(2)*u2);
% hdot2 = jacobian(h, x2)*(f22+sys.gg(2)*u2);
%% Constraint:
pconstr_12 = L1 >= 0;
% pconstr_22 = sos_constrain_2(hdot,L1,h);
pconstr_22 = hdot - L1*h >= 0;
pconstr = [pconstr_12; pconstr_22];
%% Solver parameters
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek'; % 'sedumi'
[info,dopt] = sosopt(pconstr,[x1;x2],opts);
%% Check output
if info.feas
    u_output = subs(u2,dopt);
    L1 = subs(L1,dopt);
else
    fprintf('Feasible SOS polynomial can not find.======\n');
    return;
end