function sys = prepare_data
%%
load test_obs_22-01-28-16-01.mat
dbstop if error
format long
syms x1 x2
%% import system
sys = system_formulate;
%% Chebyshev interpolants value
sz = 3;
poly_deg = 4;
it = 800;
noise = 5;
deg = 1;
idx = 1:145;
dxdt2 = 0;
X = [1];
%% The second dimension
sys.f1 = x2;
f2_appro = sos_cheb_controller(deg,sz);
%% Set up training set
train_state = data(idx,:);
train_input = action(idx,:);
%% Setup testing set
test_state = data;
test_state(idx,:)=[];
test_input = action;
test_input(idx,:) = [];
%%
x_ = [train_state(:,1) train_input(:,1)];
xtest_ = [asin(test_state(:,1)) test_input(:,1)];
%% Compute \dot_{x_2}
y_o = -sys.g/sys.l*train_state(:,1); y_bak = double(subs(f2_appro,{x1,x2},{train_state(:,1),train_input(:,1)}));
y_ = y_o - y_bak;
ytest_o = -sys.g/sys.l*test_state(:,1); ytest_bak = double(subs(f2_appro,{x1,x2},{test_state(:,1),test_input(:,1)}));
ytest_ = ytest_o - ytest_bak;
%% Gaussian Process
[mean2, ~] = gpr_xdot2(x_,y_,xtest_,ytest_,it,noise,poly_deg);
pvar x1 x2
X1 = p2s(monomials(x1,(0:poly_deg)));
X2 = p2s(monomials(x2,(0:poly_deg)));
for i = 2:length(X1)
    Y = [X; X1(i); X2(i)];
end
for i = 1:length(Y)
    dxdt2 = vpa(dxdt2 + Y(i)*mean2(i));
end
sys.f2 = dxdt2 + f2_appro;
sys.f = [p2s(x2); sys.f2];
end