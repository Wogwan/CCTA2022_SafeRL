function u_output = sos_min_u_sosp_nosim(sys, L1, f2)
% coder.extrinsic('sym');
% coder.extrinsic('syms');
% coder.extrinsic('pvar');
% coder.extrinsic('struct');
% coder.extrinsic('nvar');
% coder.extrinsic('polynomial');
% coder.extrinsic('polydecvar');
% coder.extrinsic('sosdecvar');
% coder.extrinsic('monomials');
% coder.extrinsic('subsref');
% coder.extrinsic('jacobian');
% coder.extrinsic('sosopt');
% coder.extrinsic('subs');
% coder.extrinsic('sos_compute');
% coder.extrinsic('sos_poly_add');
% coder.extrinsic('sos_fixed_barrier');
% coder.extrinsic('sos_calculate_hdot2');
% coder.extrinsic('sos_constrain_2');
%%
% x1=polynomial; x2=polynomial; 
% u2=polynomial; uc2=polynomial; 

% pconstr=polynomial; pconstr_12=polynomial; pconstr_22=polynomial; 
% info.feas=0; dopt=struct; opts.from=''; opts.solve='';
% u_degree=2;
% h=polynomial; i_nr = 0;
%%
pvar x1 x2 obj_c
f1 = x2;
%% Polynomial
[u2,uc2] = polydecvar('u_w2',monomials([x1;x2],0:2)); 
% h = sos_fixed_barrier;
% hdot = sos_calculate_hdot2(h,f1,f2,u2);
h = 1 - (x1^2+x1*x2+x2^2);
hdot = jacobian(h, x1)*(f1) + jacobian(h, x2)*(f2+u2);
%% Constraint 1:
pconstr = [];
% for i_nr = 1:length(uc2)
%     pconstr_mid = [uc2(i_nr) >= -sys.input_limit; uc2(i_nr) <= sys.input_limit];
%     pconstr = [pconstr; pconstr_mid];
% end
%%
% pconstr_22 = hdot - L1*h >= 0;
%%
pconstr_22 = hdot - L1*h - obj_c >= 0;
pconstr_33 = obj_c >= 0;
pconstr = [pconstr;pconstr_22];
% pconstr = [pconstr;pconstr_22;pconstr_33];
%% Constraint 2:
% pconstr_11 = uc2 >= -sys.input_limit;
% pconstr_12 = uc2 <= sys.input_limit;
% pconstr_22 = hdot - L1*h >= 0;
% pconstr = [pconstr_11;pconstr_12;pconstr_22];
%% Objective
% obj = sum(uc2);
% obj = trace(uc2);
% obj = obj_c;
obj = -obj_c;
%% Solver parameters
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek'; %'sedumi'
% [info,dopt] = sosopt(pconstr,[x1;x2],opts);
[info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
%% Check output
if info.feas
    u_output = subs(u2,dopt);
    obj = subs(obj_c,dopt);
else
    u_output = false;
    fprintf('Minimum input can not find.======\n');
    return;
end