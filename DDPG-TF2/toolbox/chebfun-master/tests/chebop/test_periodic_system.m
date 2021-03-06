function pass = test_periodic_system(pref)
% Test 'periodic' syntax for system of ODEs.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 2e-10;

%% System of nonlinear ODEs.
%  u - v' + v = 0, u'' - cos(v) = cos(x).

% Define the domain.
dom = [-pi pi];

% Define the rhs, and the intial guesses.
f = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];
u0 = f;

% Define the non-linear operator.
N = chebop(@(x, u, v) [ u - diff(v) + v ; diff(u, 2) - cos(v) ], dom);

% Solve using the periodic tag.
N.bc = 'periodic';
N.init = u0;
pref.discretization = @chebcolloc2;
u = solvebvp(N, f, pref);

% Solve imposing directly the periodic boundary condition.
N.bc = @(x, u, v) [ u(dom(1)) - u(dom(2)); ...
                    feval(diff(u), dom(1)) - feval(diff(u), dom(2)) ; ...
                    v(dom(1)) - v(dom(2)) ];
N.init = u0;
v = N \ f;

% Compare.
err(1) = norm(u - v, inf);

%% Test the TRIGCOLLOC class. FIRST AND SECOND ORDER LINEAR ODEs.
% u - v' = 0, u'' + v = cos(x), on [-pi pi].

% Set domain, operator L, and rhs f.
dom = [-pi, pi];
L = chebop(@(x, u, v) [ u - diff(v) ; diff(u, 2) + v ], dom);
L.bc = 'periodic';
F = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];

% Solve with TRIGCOLLOC.
U = L \ F;

u = U{1}; 
v = U{2}; 
err(2) = norm(u - diff(v));
err(3) = norm(diff(u, 2) + v - F{2});
err(4) = abs(u(dom(1)) - u(dom(2)));
err(5) = abs(feval(diff(u), dom(1)) - feval(diff(u), dom(2)));
err(6) = abs(v(dom(1)) - v(dom(2)));
err(7) = abs(feval(diff(v), dom(1)) - feval(diff(v), dom(2)));
classCheck(1) = isequal(get(u.funs{1}, 'tech'), @trigtech) && ...
     isequal(get(v.funs{1}, 'tech'), @trigtech);

%% Test the TRIGCOLLOC class. FIRST AND SECOND ORDER NONLINEAR ODEs.
%  u - v' + v = 0, u'' - cos(v) = cos(x).

% Set domain, operator L, and rhs f.
dom = [-pi, pi];
N = chebop(@(x, u, v) [ u - diff(v) + v ; diff(u, 2) - cos(v) ], dom);
N.bc = 'periodic';
F = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];
N.init = F;

% Solve with TRIGCOLLOC.
U = N \ F;

u = U{1}; 
v = U{2}; 
err(8) = norm(u - diff(v) + v);
err(9) = norm(diff(u, 2) - cos(v) - F{2});
err(10) = abs(u(dom(1)) - u(dom(2)));
err(11) = abs(feval(diff(u), dom(1)) - feval(diff(u), dom(2)));
err(12) = abs(v(dom(1)) - v(dom(2)));
err(13) = abs(feval(diff(v), dom(1)) - feval(diff(v), dom(2)));
classCheck(2) = isequal(get(u.funs{1}, 'tech'), @trigtech) && ...
     isequal(get(v.funs{1}, 'tech'), @trigtech);

%% Test the TRIGSPEC class. FIRST AND SECOND ORDER LINEAR ODEs.
% u - v' = 0, u'' + v = cos(x), on [-pi pi].

% Set domain, operator L, and rhs f.
dom = [-pi, pi];
L = chebop(@(x, u, v) [ u - diff(v) ; diff(u, 2) + v ], dom);
L.bc = 'periodic';
F = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];

% Solve with TRIGSPEC.
pref.discretization = @trigspec;
U = solvebvp(L, F, pref);

u = U{1}; 
v = U{2}; 
err(14) = norm(u - diff(v));
err(15) = norm(diff(u, 2) + v - F{2});
err(16) = abs(u(dom(1)) - u(dom(2)));
err(17) = abs(feval(diff(u), dom(1)) - feval(diff(u), dom(2)));
err(18) = abs(v(dom(1)) - v(dom(2)));
err(19) = abs(feval(diff(v), dom(1)) - feval(diff(v), dom(2)));
classCheck(3) = isequal(get(u.funs{1}, 'tech'), @trigtech) && ...
     isequal(get(v.funs{1}, 'tech'), @trigtech);


%%
pass = [classCheck, err < tol];

end
