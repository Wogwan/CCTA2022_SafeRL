function sos_python_cbf_main(a_rl, in_1, in_2)
a_rl = 0.1; input_1 = 0.1; input_2 = 0.1;
pvar x1 x2
sys.g = 10; sys.gg = [0;1];
sys.l = 1.; % Real value = 1.0
sys.m = 1.; % Real value = 1.0
sys.u_degree = 2; sys.L_au = 6;
sys.input_limit = 8;
% Taylor based equation: taylor(10*sin(x1))
f1 = x2;
f2 = x1^5/12 -5/3*x1^3+10*x1+a_rl; % f2 = 10*sin(x1)+a_rl
h = 1 - (x1^2+x1*x2+x2^2);
%
[L1,L1_Q] = sosdecvar('L1_w',monomials([x1;x2],0:sys.L_au/2));
[u2,uc2] = polydecvar('u_w2',monomials([x1;x2],0:sys.u_degree)); 
hdot = jacobian(h, x1)*(f1) + jacobian(h, x2)*(f2+u2);
%
pconstr = [];
for i_nr = 1:length(uc2)
    pconstr_mid = [uc2(i_nr) >= -sys.input_limit; uc2(i_nr) <= sys.input_limit];
    pconstr = [pconstr; pconstr_mid];
end
pconstr_11 = L1 >= 0;
pconstr_12 = hdot - L1*h >= 0;
pconstr = [pconstr;pconstr_11;pconstr_12];
%
% obj=u2;
obj = sum(uc2);
% Solver parameters
opts = sosoptions;
opts.form = 'kernel';
opts.solver = 'mosek'; %'sedumi'
[info,dopt] = sosopt(pconstr,[x1;x2],obj,opts); % [info,dopt] = sosopt(pconstr,[x1;x2],opts);
%% Check output
if info.feas
    u_output = subs(u2,dopt);
else
    u_output = false;
    fprintf('Minimum input can not find.======\n');
    return;
end

end