function sys = prepare_poly(k)
%%
dbstop if error
format long
pvar x1 x2
%% import system
sys = system_formulate;
sys.f2 = dxdt2 + f2_appro;
sys.f = [x2 
        sys.f2];
end