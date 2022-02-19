function [mean1,hyp1] = gpr_xdot2(x,y,xtest,ytest,it,noise,poly_deg)
%%
tic
syms x1 x2;
dxdt2 = 0;
X = [1];
%% Mean function
m1 = {@meanPoly,poly_deg};
hyp_m1 = ones([2*poly_deg 1]);
m2 = {@meanConst};
hyp_m2 = 0;
meanfunc = {'meanSum',{m2,m1}}; hyp.mean = [hyp_m2; hyp_m1];
%%
% sf = 0.01; ell = 0.08;
sf = 0.001; ell = 0.002;
cov1 = {@covSEiso}; hyp_cov1 = log([ell;sf/2]);
covfunc = cov1; hyp.cov = hyp_cov1;
%% Lik and Inf function
likfunc = @likGauss; sn = noise; hyp.lik = log(sn);
inf_func = {'infGaussLik'};
hyp1 = minimize(hyp, @gp, -it, inf_func, meanfunc, covfunc, likfunc, x, y);
%%
[ymu, ys2, fmu, fs2] = gp(hyp1, inf_func, meanfunc, covfunc, likfunc, x, y, xtest);
mean1 = hyp1.mean;
toc
%%
pvar x1 x2
X1 = p2s(monomials(x1,(1:poly_deg)));
X2 = p2s(monomials(x2,(1:poly_deg)));
for i = 1:length(X1)
    X = [X; X1(i); X2(i)];
end
for i = 1:length(X)
    dxdt2 = dxdt2 + X(i)*mean1(i);
end
%% Plot
% m = size(x,1); n = size(xtest,1);
% figure(801);clf;hold on;
% error = abs(ytest-ymu);
% rmse = sqrt(sum(error.^2)/n)
% plot((1:n)',ymu,'k*','LineWidth',0.8);
% plot((1:n)',ytest, 'r+', 'LineWidth',0.7);
% figure(802);clf;hold on;
% syms x1 x2
% dXtr_3 = reshape(double(subs(dxdt2,{x1,x2},{xtest(:,1),xtest(:,2)})),1,[])';
% error_2 = abs(ytest-dXtr_3);
% rmse_poly = sqrt(sum(error_2.^2)/n)
% plot((1:n)',dXtr_3,'bo','LineWidth',0.7);
% plot((1:n)',ytest,'r+','LineWidth',0.7);
end