V = polydecvar('cv',zV );
Vdot = jacobian(V,x)*f;
[sV,csV] = sosdecvar('csv', zVgram );
[s1,cs1] = sosdecvar('cs',z1);

sosc = polyconstr;
sosc(1) =  -( (Vdot+L2)+(g-V)*s2 ) >=0;
% sosc(2) =  -( (V-g)+(b-p)*s1 ) >=0 ;
% sosc(3) = s1 >= 0 ;
% sosc(4) = V-L1==sV;
% sosc(5) = sV>=0;

sosc(2) = s1 >= 0 ;
sosc(3) = V-L1==sV;
sosc(4) = sV>=0;


sopts = sosoptions;
%sopts.simplify = 'off';
[info,dopt,sossol] = sosopt(sosc,x,sopts);

subs(cs1,dopt)
sossol(2).Q

subs(csV,dopt)
sossol(4).Q

