format long
hold on;
syms x1 x2;
dom = 1.5;
g = 10;
l = 1;
m = 1;
%%
figure(22);clf;hold on;
f = @(t,x)[
    x(2)
    -10.259792327880859375*x(1)^4+0.25691540874258537533828909242099*x(1)-0.00425688736140727996826171875*x(2)^4+0.00032671549706708713500802332418971
    ];
%     x(2)
%     -0.004349947907030582427978515625*x(1)^4-0.005278733558952808380126953125*x(1)^3+0.000331557705067098140716552734375*x(1)^2-9.9989249945429218985326613733378*x(1)+0.08957691490650177001953125*x(2)^4-0.089119434356689453125*x(2)^3+0.01762610487639904022216796875*x(2)^2+0.00110464193858206272125244140625*x(2)-0.00014116639795298813186974484779057
%     ];
% f = @(t,x)[
%     x(2)
%      -g/l*sin(x(1))
%     ];

%% Field [-2 2]
% vectfield(f,-dom:dom/10:dom,-dom:dom/10:dom); hold on; % field [-5 5 -5 5]
for xd = -dom:dom/5:dom
    for yd = -dom:dom/5:dom
        [ts,ys] = ode45(f,[0,20],[xd;yd]);
        h1 = plot(ys(:,1),ys(:,2), 'k');hold on;
        refreshdata; drawnow;
    end
end
xlim([-dom dom]); ylim([-dom dom]);
plot(0,0,'c*','linewidth',2);hold on;
set(gca, 'LooseInset', [0,0,0,0]);
title('');