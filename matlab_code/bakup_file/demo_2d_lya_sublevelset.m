clear;
tic
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
f = [x2-x1
    -0.39306944559376340297962570957679*x1^4+0.30314332186361650500749931325117*x1^3+x1^2*x2+0.46065476097217313011800143840446*x1^2-0.45469239127041680137431001185178*x1+0.049252630679030955096475707932768*x2^4+0.045789722029436145944725211620607*x2^3-0.063117281796006285965461302112089*x2^2+0.009401592216201428570121478855981*x2+0.0023141262623443403789735839382047
    ];
gg = [1;1];
input = [gg(1)*u1;gg(2)*u2];
% V = 1*x1^4+2*x2^4+2*x1^2*x2^2+1*x1^2+1*x2^2+1*x1*x2;
V = 1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1;
%%
C0 = 0.1;
cc = 1;
k = 4;
k_l = 4;
dom = 15;
domain = [-dom dom -dom dom];
%%
C1 = (x1+4)^2+(x2-5)^2-4;
C2 = (x1+0)^2+(x2+5)^2-4;
C3 = (x1-5)^2+(x2-0)^2-4;
%%
C = [C1;C2;C3];
kk = 1;
%%
figure_id = 11;
figure(figure_id);clf;hold on;
[~,~]=pcontour(C(1),0,domain,'r');
[~,~]=pcontour(C(2),0,domain,'r');
[~,~]=pcontour(V,C0,domain,'g');
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'r');
end
%%
solU = [];
v_c = [];
iter = 1;
%%
while double(cc)-double(C0) >= 1e-6
    iter = iter + 1;
    if iter ~= 1
        C0 = cc;
    end
    [solu1,solu2,solL,kk] = sos_function_v(f,gg,k,k_l,V,C0);
    if kk == 0
        break
    end
    [cc,kk,solu1,solu2] = sos_function_v2(f,gg,k,k_l,V,C,dom,solL,figure_id);
    v_c = [v_c; double(cc)];
    solU = [solU;[solu1,solu2]];
    if kk == 0
        figure(figure_id);hold on;
        [~,~]=pcontour(V,v_c(end),domain,'b');
        break
    end
end
figure(figure_id);hold on;
[~,~]=pcontour(V,v_c(end),domain,'b');
%%
fprintf('Sublevel set of V0(x) is %.14f. \n',v_c(end));
toc