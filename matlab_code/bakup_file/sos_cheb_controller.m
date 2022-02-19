function f = sos_cheb_controller(deg,sz)
%%
sz = 2;
deg = 5;
a = -sz; b = sz;
y = chebfun('10*sin(x)',[a,b],'splitting','on'); % Modified here with different non-polynomial g
y_3 = minimax(y,deg); c_3 = chebcoeffs(y_3);
% Plot the chebyshev approximation result with original function
syms x1;
T = chebyshevT([0:deg],x1);
if length(T)~= length(c_3)
    c_3 = [c_3; 0];
    f0 = vpa(T*c_3);
elseif length(T)==length(c_3)
    f0 = vpa(T*c_3);
end
x2 = x1/sz;
f = subs(f0,x1,x2);
figure(3);clf;
plot(y,'r-'); xlim([a b]);ylim([-10 10]);hold on;
ezplot(f); hold on;
plot(y_3,'g--');

%% Taylor
end

