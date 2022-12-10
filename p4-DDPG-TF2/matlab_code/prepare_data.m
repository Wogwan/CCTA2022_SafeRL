function sys = prepare_data(k)
% k = 1
%%
dbstop if error
format long
syms x1 x2
%% Chebyshev interpolants value
sz = 2;
poly_deg = 3;
it = 600;
noise = 1;
deg = 3;
%%
if k == 1
    load obs_env_record.mat
    data = data_base;
    action = action_base;
    f2_appro = sos_cheb_controller(deg,sz);
else
    %% data, action, data_base, action_base, rmse_poly_cur, mean2_cur, rmse_poly_base, mean2_base
    load obs_env.mat
    load line_search_GP.mat
    load gp_model.mat
end


%% If data is not good enough, employee the previous data
% data_base = data;
% action_base = action;
% save('C:\ASCC_2022_SafeRL\utest\ASCC2022_SafeRL\DDPG-TF2\res\mat\obs_env_record.mat','action_base','data_base')

%% import system
sys = system_formulate;
length_size = length(data);
idx = 1:length_size/5*4;
dxdt2 = 0;
f2_cur = f2_appro;
X = [1];
%% The second dimension
sys.f1 = x2;
%% Set up training set
train_state = data(idx,:);
train_input = action(idx,:);
%% Setup testing set
test_state = data;
test_state(idx,:) = [];
test_input = action;
test_input(idx,:) = [];
%%
x_ = train_state(:,1);
xtest_ = test_state(:,1);
%%
if k==1
    %% Compute \dot_{x_2}
    y_o = -10*sin(train_state(:,1));
    y_bak = double(subs(f2_appro,{x1},{train_state(:,1)}));
    y_ = y_o - y_bak;
    ytest_o = -10*sin(test_state(:,1));
    ytest_bak = double(subs(f2_appro,{x1},{test_state(:,1)}));
    ytest_ = ytest_o - ytest_bak;
else
    y_o = -10*sin(train_state(:,1));
    y_bak = double(subs(f2_appro,{x1,x2},{train_state(:,1),train_state(:,2)}));
    y_ = y_o - y_bak;
    ytest_o = -10*sin(test_state(:,1));
    ytest_bak = double(subs(f2_appro,{x1,x2},{test_state(:,1),test_state(:,2)}));
    ytest_ = ytest_o - ytest_bak;
end
%% Gaussian Process
[mean2_mid, ~, rmse_poly_mid] = gpr_xdot2(x_,y_,xtest_,ytest_,it,noise,poly_deg);
%%
if k == 1
    rmse_poly_cur = rmse_poly_mid;
    mean2_cur = mean2_mid;
    mean2 = mean2_mid;
    save('..\res\mat\line_search_GP.mat','rmse_poly_cur','mean2_cur')
    sys.stop = 0;
else
    if rmse_poly_mid <= rmse_poly_cur
        rmse_poly_cur = rmse_poly_mid;
        mean2_cur = mean2_mid;
        mean2 = mean2_mid;
        save('..\res\mat\line_search_GP.mat','rmse_poly_cur','mean2_cur')
        sys.stop = 0;
    else
        mean2 = mean2_cur;
        sys.stop = 1;
    end
end
%%
pvar x1 x2
X1 = p2s(monomials(x1,(0:poly_deg)));
X2 = p2s(monomials(x2,(0:poly_deg)));
for i = 2:length(X1)
    % Y = [X; X1(i); X2(i)];
    Y = [X; X1(i)];
end
for i = 1:length(Y)
    dxdt2 = vpa(dxdt2 + Y(i)*mean2(i));
end
sys.f2 = dxdt2 + f2_appro;
sys.f = [p2s(x2); sys.f2];
f2_appro = sys.f2;
save('..\res\mat\gp_model.mat','f2_appro')

end