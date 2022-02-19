% function sys = s3_opt_lyapunov(sys)
pvar x1 x2 u1 u2 htol epsi;
format long
x = [x1;x2];
%%
dom = 10; domain = [-dom dom -dom dom];
f = [x2-x1
    -0.39306944559376340297962570957679*x1^4+0.30314332186361650500749931325117*x1^3+x1^2*x2+0.46065476097217313011800143840446*x1^2-0.45469239127041680137431001185178*x1+0.049252630679030955096475707932768*x2^4+0.045789722029436145944725211620607*x2^3-0.063117281796006285965461302112089*x2^2+0.009401592216201428570121478855981*x2+0.0023141262623443403789735839382047
    ];
gg = [1; 1];
sys = [f(1)+gg(1)*u1; f(2)+gg(2)*u2];
C1 = (x1+4)^2+(x2-5)^2-4;
C2 = (x1+0)^2+(x2+5)^2-4;
C3 = (x1-5)^2+(x2-0)^2-4;
C = [C1;C2;C3];
V0 = 1*x1^4+1*x1^3*x2+1*x1^2*x2^2+1*x1*x2^3+1*x2^4+1*x1^3+1*x1^2*x2+1*x1*x2^2+1*x2^3+1*x1^2+1*x1*x2+1*x2^2+1*x1+1*x2+1;
CC0 = 57.23191152365694;
u1 = 112.37400271765102388599188998342*x2 - 50.040555105546509651048836531118*x1 - 0.0000014928903836983044263534954956829*x1^2*x2^2 - 76.733210219968043475091690197587*x1*x2 - 772.48097300437927970051532611251*x1*x2^2 - 162.07557471662337889029004145414*x1^2*x2 - 0.0000084441309383659872303887206079764*x1*x2^3 + 0.0000081933097066728891434781584246494*x1^3*x2 + 38.658043521953828758341842330992*x1^2 - 43.911527505085800271444895770401*x1^3 - 220.01710526604153983498690649867*x2^2 - 0.00000053134452790163867559120524849958*x1^4 - 295.53486232374842757053556852043*x2^3 + 0.0000077452709307552725862943548973405*x2^4 - 29.868425418772179824600243591703;
u2 = 116.85988100511805498626927146688*x1 - 66.06755535749654484334314474836*x2 - 0.0000050288140412395193555561943854482*x1^2*x2^2 - 194.62574107076193286047782748938*x1*x2 - 267.78834733834349890457815490663*x1*x2^2 - 732.55532954320972294226521626115*x1^2*x2 + 0.0000061282118915948484533478927971384*x1*x2^3 - 0.0000037935843498656611732632731870396*x1^3*x2 - 90.067549213125758456044422928244*x1^2 - 192.33960066656609910751285497099*x1^3 + 9.0341654196018055245076538994908*x2^2 + 0.39307852943896876007912055683846*x1^4 - 52.394704864320068793404061580077*x2^3 - 0.049252424236205412377831436288034*x2^4 - 29.304083970311982199064004817046;
B = 205.39988836883230760577134788036*x1^2*x2^2 - 15.733076483919825605539699608926*x2 - 32.343368076002114719358360162005*x1 - 16.702690567188692938316307845525*x1*x2 + 58.06348388736881105387510615401*x1*x2^2 - 0.88882356086419722629443640471436*x1^2*x2 + 27.919011189967342545514839002863*x1*x2^3 + 2.3432143317831926054850555374287*x1^3*x2 - 16.00443066513419054786027118098*x1^2 - 4.0245318441951161148040227999445*x1^3 + 21.787467874077073304306395584717*x2^2 - 119.2190088341546641004242701456*x1^4 + 0.47139802020367055357397134685016*x2^3 - 139.72095091473158845474245026708*x2^4 + 10005.142017583266351721249520779;
%%

%%
c0 = 1; cc = 1.1; epsi = 1e-6; figure_id = 111;
%%
figure(figure_id);clf;hold on;
[~,~]=pcontour(C(1),0,domain,'r');
[~,~]=pcontour(C(2),0,domain,'r');
if length(C) == 3
    [~,~]=pcontour(C(3),0,domain,'r');
elseif length(C) == 4
    [~,~]=pcontour(C(3),0,domain,'r');
    [~,~]=pcontour(C(4),0,domain,'r');
else
    fprintf('The constraint number does not match.======\n');
end
hold on; [~,~]=pcontour(B,0,domain,'g');
[~,~]=pcontour(V0,double(CC0),domain,'m');
%% Hyperparameters of the SOSP @ CBF -> V
V_us = 4; V_au = 4; V_degree = 4; gamma = 0; k_u_V = 2; k_l_au = 4;
%%
kk = 1; OO = 0;
dom_2 = 10000; domain_2 = [-dom_2 dom_2 -dom_2 dom_2]; N_Lya = [];
%%
[V, kk] = sos_optimal_V1(f,gg,B,u1,u2,V_au,V_us,V_degree,C,gamma);
%%
if kk == 0
    fprintf('Suitable Lyapunov function can not find.======\n');
end
N_Lya = [N_Lya;V];
%% TEST FOR Sublevel Set
[a1,b1] = coeffs(p2s(V));
ccc = double(vpa(a1(end)));
C0 = ccc; cc = ccc+1;
%%
figure(figure_id);hold on;
[~,~]=pcontour(V,double(C0),domain,'r');
solU = []; v_c = []; iter = 0;
%%
while double(cc)-double(C0) >= epsi
    iter = iter + 1;
    if iter ~= 0
        C0 = cc;
    end
    [solL,kk]= sos_optimal_v2_origin(f,gg,k_u_V,k_l_au,V,cc);
    if kk == 0
        break
    end
    [cc,kk,solu1,solu2] = sos_optimal_v3_origin(f,gg,k_u_V,k_l_au,V,C,dom,solL,ccc,figure_id);
    v_c = [v_c; double(cc)];
    solU = [solU;[solu1,solu2]];
    if kk == 0
        figure(figure_id);hold on;
        [~,~]=pcontour(V,v_c(end),domain,'b');
        break
    end
end
%% Start to compute the control barrier function
c_b = v_c(end);
sol_B = c_b - N_Lya(end);
hold on;[~,~]=pcontour(sol_B,0,domain,'k');
V = N_Lya(end);
% end
%%
% fprintf('Sublevel set of V∗(x) is \n%s \n\n',char(vpa(p2s(V))));
% fprintf('Sublevel set of V∗(x) is %.14f. \n',v_c(end));