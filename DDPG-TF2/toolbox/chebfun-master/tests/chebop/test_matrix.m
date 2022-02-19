function pass = test_matrix(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

tol = 1e-10;

%%
D = chebop(@(u) diff(u));

D13 = [ ...
  -2.488033871712585   2.821367205045918  -0.511966128287415   0.178632794954082
   0.166666666666667  -1.333333333333333   1.333333333333333  -0.166666666666667
  -0.178632794954082   0.511966128287415  -2.821367205045918   2.488033871712584
     ];
err(1) = norm(matrix(D, 3) - D13);

D4old = [ ...
  -3.166666666666667   4.000000000000000  -1.333333333333333   0.500000000000000
  -1.000000000000000   0.333333333333333   1.000000000000000  -0.333333333333333
   0.333333333333333  -1.000000000000000  -0.333333333333333   1.000000000000000
  -0.500000000000000   1.333333333333333  -4.000000000000000   3.166666666666667
     ];
err(2) = norm(matrix(D, 4, 'oldschool') - D4old);

%%
D2 = chebop(@(u) diff(u, 2));

D22 = [ ...
   4.161760458079525  -6.990187582825714   4.323520916159047  -1.495093791412856
  -1.495093791412856   4.323520916159048  -6.990187582825714   4.161760458079524
     ];
err(3) = norm(matrix(D2, 2) - D22);

D24old = [ ...
   5.333333333333335  -9.333333333333336   6.666666666666668  -2.666666666666667
   3.333333333333334  -5.333333333333334   2.666666666666667  -0.666666666666667
  -0.666666666666666   2.666666666666667  -5.333333333333334   3.333333333333334
  -2.666666666666667   6.666666666666668  -9.333333333333336   5.333333333333335
     ];
err(4) = norm(matrix(D2, 4, 'oldschool') - D24old);

%% Boundary conditions:
D2 = chebop(@(u) diff(u, 2), [-1 1], 0);
D2bc = [1 0 0 0 ; 0 0 0 1 ; D22];
err(5) = norm(matrix(D2, 2) - D2bc);

D24oldbc = [1 0 0 0 ; D24old(2:3,:) ; 0 0 0 1]; 
err(6) = norm(matrix(D2, 4, 'oldschool') - D24oldbc);

%% Continuity conditions:
D2 = chebop(@(u) diff(u, 2), [-1 0 1]);
D22cont = [0 0 0 1 -1 0 0 0 ;
           -1 8/3 -8 19/3 19/3 -8 8/3 -1 ;
           4*D22 0*D22 ; 
           0*D22 4*D22];
err(7) = norm(matrix(D2, [2 2]) - D22cont);

D24oldcont = [0 0 0 1 -1 0 0 0 ;
              4*D24old(2:4,:) zeros(3,4) ; 
              zeros(3,4) 4*D24old(1:3,:) ; 
              -1 8/3 -8 19/3 19/3 -8 8/3 -1];
err(8) = norm(matrix(D2, [4 4], 'oldschool') - D24oldcont);        
              
%%

pass = err < tol;
 
end
