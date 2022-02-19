close all;clear;
% sosaddpath
syms x1 x2;
%% Collect data
sys = prepare_data;
%% Check Lyapunov sublevel set function
sys = s1_lya_sublevelset(sys);
%% Check permissive barrier function and optimal output
sys = s2_opt_barrier(sys);
%% Check permissive Lyapunov function
sys = s3_opt_lyapunov(sys);
