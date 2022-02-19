u = 0.5;
input_1 = 0.8;
input_2 = 3;
sys = sys_launch(u);
%%
while 1
    sys.input_limit = sys.input_limit - 0.1;
    sys.input_limit
    %% Find optimal L1
    L1 = sos_feasu(sys);
    %% Find minimize U
    u_out = sos_min_u_sosp(sys, L1);
    try
        if class(u_out) == 'polynomial'
            u_store = u_out;
            continue
        end
    catch
        u_store
        sys.input_limit
        break
    end
end
syms x1 x2
u_sos = p2s(u_store);
y = subs(u_sos, {x1,x2}, {input_1,input_2});