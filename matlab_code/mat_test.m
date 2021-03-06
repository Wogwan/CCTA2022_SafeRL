%%
dbstop if error
format long
syms x1 x2
%% Chebyshev interpolants value
sz = 2;
poly_deg = 3;
it = 600;
noise = 10;
deg = 3;
%%
load C:\ASCC_2022_SafeRL\utest\ASCC2022_SafeRL\DDPG-TF2\res\mat\obs_env_record.mat
data = data_base;
action = action_base;
f2_appro = sos_cheb_controller(deg,sz);
%% import system
sys = system_formulate;
length_size = length(data);
idx = 1:length_size/5*4;
dxdt2 = 0;
X = [1];
%% The second dimension
sys.f1 = x2;
%% Set up training set
train_state = data(idx,:);
train_input = action(idx,:);
%% Setup testing set
test_state = data;
test_state(idx,:)=[];
test_input = action;
test_input(idx,:) = [];
%%
x_ = train_state;
xtest_ = test_state;
%%
%% Compute \dot_{x_2}
y_o = sin(train_state(:,1));
y_bak = double(subs(f2_appro,{x1},{train_state(:,1)}));
y_ = y_o - y_bak;
ytest_o = test_state(:,1);
ytest_bak = double(subs(f2_appro,{x1},{test_state(:,1)}));
ytest_ = ytest_o - ytest_bak;
%% Gaussian Process
[mean2_mid, ~, rmse_poly_mid] = gpr_xdot2(x_,y_,xtest_,ytest_,it,noise,poly_deg);
%%
rmse_poly_cur = rmse_poly_mid;
mean2_cur = mean2_mid;
mean2 = mean2_mid;
save('C:\ASCC_2022_SafeRL\utest\ASCC2022_SafeRL\DDPG-TF2\res\mat\line_search_GP.mat','rmse_poly_cur','mean2_cur')
sys.stop = 0;

%%
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
f2_appro = sys.f2;
%%
% save('C:\ASCC_2022_SafeRL\utest\ASCC2022_SafeRL\DDPG-TF2\res\mat\gp_model.mat','f2_appro')

% end