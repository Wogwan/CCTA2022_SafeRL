clear
%% System information
sys.dom = 6; 
sys.domain = [-sys.dom sys.dom -sys.dom sys.dom];
sys.dt = 0.05;
sys.g = 10;
sys.l = 1.4; % Real value = 1.0
sys.m = 1.4; % Real value = 1.0
sys.h_degree = 2;

%% Build ode
syms x1 x2
% Input current state [theta_t, theta_dot_t]
% x1_limit = 1;
% x1_t = -x1_limit + 2*x1_limit*rand(1,1);
% x2_limit = 8;
% x2_t = -x2_limit + 2*x2_limit*rand(1,1);
% x_t = [x1_t; x2_t];

% state: -0.9259530516281068 -1.2141997278723683e-06
% action_rl: [1.78829904] |  u_BAR: [-0.49779296] | u_bar: [2.70548417] | output: [3.99599025]
x1_t = -0.9259530516281068;
x2_t = -1.2141997278723683e-06;
x_t = [x1_t; x2_t];
u_rl = 1.78829904;

% G term
G = [3/(sys.m*sys.l^2)*sys.dt^2
    3/(sys.m*sys.l^2)*sys.dt];

% Prepare the t+1 state value [theta_{t+1}, theta_dot_{t+1}]
x1_t_ = -3*sys.g/(2*sys.l)*sin(x1_t+pi)*sys.dt^2 + x2_t*sys.dt + x1_t + 3/(sys.m*sys.l^2)*u_rl*sys.dt^2; 
x2_t_ = x2_t - 3*sys.g/(2*sys.l)*sin(x1_t+pi)*sys.dt + 3/(sys.m*sys.l^2)*u_rl*sys.dt;
x_t_ = [x1_t_; x2_t_];

%%
% min a_t
% s.t. h(s_{t+1})-h(s_{t}) >= 0
%      h >= 0

% h
% a_t

% input value
x_t_input = [1; x1_t; x2_t; x1_t^2; x1_t*x2_t; x2_t^2];
x_t__input = [1; x1_t_; x2_t_; x1_t_^2; x1_t_*x2_t_; x2_t_^2];
%% Build SOS minimize input : compute a minimize control based on a given barrier function 
pvar x1 x2 u_bar
x = [x1;x2];
[h1,h1_Q] = polydecvar('h1_w',monomials(x,0:sys.h_degree)); 
% [h1,h1_Q] = sosdecvar('h1_w',monomials(x,0:sys.h_degree/2)); 
% Constraint:
pconstr_1 = sum(x_t__input.*h1_Q) + G(2)*u_bar - sum(x_t_input.*h1_Q) >= 0;
pconstr_2 = h1 >= 0;
pconstr = [pconstr_1; pconstr_2];
% Objective
obj = -u_bar;
% Solver parameters
opts = sosoptions;
opts.form = 'kernel';
% opts.solver = 'mosek';
opts.solver = 'sedumi';
[info,dopt] = sosopt(pconstr,[x1;x2],obj,opts);
% Check output
if info.feas
    u_output = subs(u_bar,dopt);
    h_output = subs(h1,dopt);
else
    fprintf('Minimum input can not find.======\n');
    return;
end