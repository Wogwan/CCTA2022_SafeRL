clear;tic;
pvar x1 x2 u htol epsi;
format long
x = [x1;x2];
dom = 10;
%%
f = [x2-x1
    -0.39306944559376340297962570957679*x1^4+0.30314332186361650500749931325117*x1^3+x1^2*x2+0.46065476097217313011800143840446*x1^2-0.45469239127041680137431001185178*x1+0.049252630679030955096475707932768*x2^4+0.045789722029436145944725211620607*x2^3-0.063117281796006285965461302112089*x2^2+0.009401592216201428570121478855981*x2+0.0023141262623443403789735839382047
    ];
gg = [1;1];
%%
% V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2; 
V = 1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1;
C0 = 57.23191152365694;
%%
C1 = (x1+4)^2+(x2-5)^2-4;
C2 = (x1+0)^2+(x2+5)^2-4;
C3 = (x1-5)^2+(x2-0)^2-4;
%%
C = [C1;C2;C3];
trace_Q1 = 1; trace_Q = 0; mm = 0; kk = 1; iter = 0;
k = ['r','g','b','m','c','k','y'];
%%
sol_B = C0 - V;
solh = sol_B;
%%
k_u = 4; k_h = 4; L_us = 4; gamma = 0;
L_au = 2; 
%L_au = 6; 
%%
figure_id = 12;
figure(figure_id+1);clf;hold on;
figure(figure_id+2);clf;hold on;
figure(figure_id);clf;hold on;
domain = [-dom dom -dom dom];
xlim([-dom dom]); ylim([-dom dom]); hold on;
[~,~]=pcontour(V,C0,domain,'b'); hold on;
[~,~]=pcontour(C1,0,domain,'k'); hold on;
[~,~]=pcontour(C2,0,domain,'k'); hold on;
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'k');         
end
axis(domain); TRACE = [];
Barrier = []; Control = []; Barrier_plus = [];
%%
while 1
    iter = iter+1
    record_Q = trace_Q
    %%
    [SOLu1,SOLu2,SOL1,SOL2,kk] = sos_function_1(f,k_u,L_au,solh,V,gamma,gg);
    if kk == 0
        break
    end
    Control = [Control; [SOLu1 SOLu2]];
    %%
    [solh,trace_Q,kk] = sos_function_2(iter,f,k_h,SOLu1,SOLu2,SOL1,SOL2,gamma,V,C,dom,gg,L_us,figure_id);
    if kk == 0
        break
    end
    TRACE = [TRACE; double(trace_Q)];
    Barrier = [Barrier; solh];
    %% Optimal the set
    kk = 1; OO = 0;
    %%
    if kk == 0
        fprintf('Advanced Barrier Function can not find.======\n');
    end
    figure(figure_id+1);clf;hold on;
    [~,~]=pcontour(V,C0,domain,'b'); hold on;
    [~,~]=pcontour(C1,0,domain,'k'); hold on;
    [~,~]=pcontour(C2,0,domain,'k'); hold on;
    [~,~]=pcontour(C3,0,domain,'k');
    [~,~]=pcontour(V,C0,domain,'r'); hold on;
    if mod(iter,7) == 0
        [~,~]=pcontour(solh,0,domain,k(7)); hold on;             % Plot the original Lyapunov sublevel set
    else
        [~,~]=pcontour(solh,0,domain,k(mod(iter,7))); hold on;             % Plot the original Lyapunov sublevel set
    end
    refreshdata; drawnow;
end
toc
A = [];
for iter = 1:length(Barrier)
    A = [A; [Control(iter,:) Barrier(iter)]];
end
%%
fprintf('Permissive B(x) is \n%s \n\n',char(vpa(p2s(A(end,3)))));
fprintf('Control Input u1(x) is \n%s \n\n',char(vpa(p2s(A(end,1)))));
fprintf('Control Input u2(x) is \n%s \n\n',char(vpa(p2s(A(end,2)))));