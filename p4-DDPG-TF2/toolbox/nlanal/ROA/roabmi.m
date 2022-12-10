function [bopt,gopt,s1opt,s2opt,Vopt] = roabmi(f,x,p,zV,opts,ic)
% function [bopt,gopt,s1opt,s2opt,Vopt] = roabmi(f,x,p,zV,opts,ic)
%
% PJS 7/2/2009   Initial coding

%------------------------------------------------------------------
% Error Checking
%------------------------------------------------------------------

%------------------------------------------------------------------
% Get options or set defaults for unspecified options
%------------------------------------------------------------------
if nargin<5
    opts = [];
    ic = [];
elseif nargin<6
    ic = [];
end

% Set default options
if isempty(opts)
    opts.L1 = [];
    opts.L2 = [];
    opts.z1 = [];
    opts.z2 = [];
    opts.solveropts = [];
end

if ~isfield(opts,'L1') || isempty(opts.L1)
    L1 = 1e-6*(x'*x);
else
    L1 = opts.L1;
end
if ~isfield(opts,'L2') || isempty(opts.L2)
    L2 = 1e-6*(x'*x);
else
    L2 = opts.L2;
end
if ~isfield(opts,'z1') || isempty(opts.z1)
    Vdeg = zV.maxdeg;
    z1maxd = ceil((Vdeg-p.maxdeg)/2);
    z1 = monomials(x, 0:z1maxd ); 
else
    z1 = opts.z1;
end
if ~isfield(opts,'z2') || isempty(opts.z2)
    z2 = monomials(x, 1:2 ); 
else
    z2 = opts.z2;
end

%------------------------------------------------------------------
% Set up ROA BMI
%------------------------------------------------------------------

% Create decision vars for level sets
pvar gam bet;

% Create decision vars for Lyapunov function and SOS multipliers
[s1,c1mat] = sosdecvar('c',z1);
[s2,c2mat] = sosdecvar('d',z2);
[V,cV] = polydecvar('e',zV);

% Stack all ROA decision variables (BMI will also include null space 
% decision vars associated with the Gram matrix constraint)
tmpu=triu( ones(size(c1mat)) ); 
idxu1 = find(tmpu); 
c1 = c1mat(idxu1);
lc1 = length(c1);
lz1 = length(z1);

tmpu=triu( ones(size(c2mat)) ); 
idxu2 = find(tmpu); 
c2 = c2mat(idxu2);
lc2 = length(c2);
lz2 = length(z2);

lcV = length(cV);
lzV = length(zV);

decvars = char([bet; gam; c1; c2; cV]);
ldv = length(decvars);

% Formulate BMI constraints:
%   A0 + sum_i xi*Ai + sum_i sum_j x_i x_j Kij <= 0
% See bmi2pen for interface.
% x = [decvars; n3; n4; n5] where n3, n4, n5 are the  null space
% variables associated with constraints 3, 4, and 5.

% ------- Constraint 1: s1 in SOS
msizes = lz1;
A0 = sparse(lz1^2,1);
A = sparse(lz1^2,ldv);
c1idx = 2+(1:lc1);
A(:,c1idx) = sparse(idxu1,1:lc1,-ones(lc1,1),lz1^2,lc1);

% ------- Constraint 2: s2 in SOS
msizes(2) = lz2;
A0 = [A0; sparse(lz2^2,1)];
Atmp = sparse(lz2^2,ldv);
c2idx = 2+lc1+(1:lc2);
Atmp(:,c2idx) = sparse(idxu2,1:lc2,-ones(lc2,1),lz2^2,lc2);
A = [A; Atmp];

% ------- Constraint 3: V-L1  in SOS
s = V-L1;
[g0,g,d] = collect(s,x);
[isdec,didx]=ismember(char(d),decvars);   
[w3,Q3,N3,D3] = gramsol(g0,g);  
Q3 = Q3(:);
lw3 = length(w3);

msizes(3) = lw3;
A0 = [A0; -Q3];
A = blkdiag(A,-N3);
A(end-lw3^2+1:end,didx) = -D3;

% ------- Constraint 4: -( (V-gamma)+(beta-p)*s1 ) in SOS 
s = -( (V-gam)+(bet-p)*s1 );
[g0,g,d] = collect(s,x);

% Separate linear / bilinear terms
lidx = find(sum(d.degmat,2)==1);
nl = length(lidx);
bidx = find(sum(d.degmat,2)==2);
nb = length(bidx);
g = [g(lidx)' g(bidx)'];
dl = d(lidx);
db = d(bidx);

% Create linear terms
[isdec,dlidx]=ismember(char(dl),decvars);   
[w4,Q4,N4,D4] = gramsol(g0,g);  
Q4 = Q4(:);
lw4 = length(w4);

msizes(4) = lw4;
A0 = [A0; -Q4];
A = blkdiag(A,-N4);
A(end-lw4^2+1:end,dlidx) = -D4(:,1:nl);

% Split bilinear variables and store bilinear terms
% (Matrix for bilinear terms will be created after we know
%  how many null space variables are needed)
db1 = cell(nb,1);
db2 = db1;
deg = db.coefficient'*db.degmat;
var = db.varname;
[ridx,cidx]=find(deg==2);
db1(ridx) = var(cidx);
db2(ridx) = var(cidx);
[cidx,ridx]=find(deg'==1);
ridx = reshape(ridx,2,length(ridx)/2);
cidx = reshape(cidx,2,length(cidx)/2);
db1(ridx(1,:)) = var(cidx(1,:));
db2(ridx(1,:)) = var(cidx(2,:));

[isdec,db1idx]=ismember(db1,decvars);   
[isdec,db2idx]=ismember(db2,decvars);   
storeK4 = {db1idx,db2idx,-D4(:,nl+1:end)};
% dbidx = sub2ind([ldv ldv],db1idx,db2idx);
% K = sparse(size(At,1),ldv^2);
% K(end-lw4^2+1:end,dbidx) =  -D4(:,nl+1:end);

% ------- Constraint 5: -( (grad(V)*f+L2)+(gamma-V)*s2 ) in SOS
Vdot = jacobian(V,x)*f;
s = -( (Vdot+L2)+(gam-V)*s2 );
[g0,g,d] = collect(s,x);

% Separate linear / bilinear terms
lidx = find(sum(d.degmat,2)==1);
nl = length(lidx);
bidx = find(sum(d.degmat,2)==2);
nb = length(bidx);
g = [g(lidx)' g(bidx)'];
dl = d(lidx);
db = d(bidx);

% Create linear terms
[isdec,dlidx]=ismember(char(dl),decvars);   
[w5,Q5,N5,D5] = gramsol(g0,g);  
Q5 = Q5(:);
lw5 = length(w5);

msizes(5) = lw5;
A0 = [A0; -Q5];
A = blkdiag(A,-N5);
A(end-lw5^2+1:end,dlidx) = -D5(:,1:nl);

% Split bilinear variables and store bilinear terms
% (Matrix for bilinear terms will be created after we know
%  how many null space variables are needed)
db1 = cell(nb,1);
db2 = db1;
deg = db.coefficient'*db.degmat;
var = db.varname;
[ridx,cidx]=find(deg==2);
db1(ridx) = var(cidx);
db2(ridx) = var(cidx);
[cidx,ridx]=find(deg'==1);
ridx = reshape(ridx,2,length(ridx)/2);
cidx = reshape(cidx,2,length(cidx)/2);
db1(ridx(1,:)) = var(cidx(1,:));
db2(ridx(1,:)) = var(cidx(2,:));

[isdec,db1idx]=ismember(db1,decvars);   
[isdec,db2idx]=ismember(db2,decvars);   
storeK5 = {db1idx,db2idx,-D5(:,nl+1:end)};
% dbidx = sub2ind([ldv ldv],db1idx,db2idx);
% K( (end+1):(end+lw5^2),dbidx) =  -D5(:,nl+1:end);

% ------- Create BMI matrix
nvars = size(A,2);
m = size(A,1);
K = sparse(m,nvars^2);
tmpidx = [0 cumsum(msizes.^2)];

ridx4 = tmpidx(4)+1:tmpidx(5);
cidx4 = sub2ind([nvars nvars],storeK4{1},storeK4{2});
K(ridx4,cidx4) = storeK4{3};

ridx5 = tmpidx(5)+1:tmpidx(6);
cidx5 = sub2ind([nvars nvars],storeK5{1},storeK5{2});
K(ridx5,cidx5) = storeK5{3};

% ------- Objective function: min (-beta)
fobj = sparse(nvars,1);
fobj(1) = -1;

% ------- Create BMI initial condition

if ~isempty(ic)
    bet0 = ic.beta;
    gam0 = ic.gamma;
    s10 = ic.s1;
    s20 = ic.s2;
    V0 = ic.V;
        
    % Get intial values for s1
    if isequal(z1,1)
        c10 = double(s10);
        feas1 = 1;
    else
        [AA,Rbasis] = buildLmat(z1);  
        bb = poly2basis(s10,Rbasis);
        KK.s = msizes(1);
        pars.fid = 0;
        [xx1,yy1,info1]=sedumi(AA,bb,[],KK,pars);
        feas1 = info1.pinf==0 && info1.dinf==0 && ...
            (info1.numerr==0 || info1.numerr==1);
        c10 = reshape( xx1 , msizes(1)*[1 1]);
        c10 = c10(idxu1);
    end
    %s10-subs(s1,c1,c10)

    % Get intial values for s2
    if isequal(z1,1)
        c20 = double(s20);
        feas2 = 1;
    else
        [AA,Rbasis] = buildLmat(z2);  
        bb = poly2basis(s20,Rbasis);
        KK.s = msizes(2);
        pars.fid = 0;
        [xx2,yy2,info2]=sedumi(AA,bb,[],KK,pars);
        feas2 = info2.pinf==0 && info2.dinf==0 && ...
            (info2.numerr==0 || info2.numerr==1);
        c20 = reshape( xx2 , msizes(2)*[1 1]);
        c20 = c20(idxu2);
    end
    %s20-subs(s2,c2,c20)

    % Get initial values for V
    cV0 = full(poly2basis(V0,zV));
    %V0-subs(V,cV,cV0)

    x0  = [bet0; gam0; c10(:); c20(:); cV0(:)];    
    x0  = [x0; zeros(nvars-ldv,1)];
    x0x0 = x0*x0';
    x0x0 = x0x0(:);

    % Solve for initial null space vars  
    % (Note that the constant term in Sedumi has opposite sign
    %  convention than that used in bmi2pen)
    cc = -(A0 + A*x0+ K*x0x0);    
    AA = A(:,ldv+1:end); 

    KK.s = msizes(3:end);
    ptr = msizes(1)^2+msizes(2)^2+1;
    cc = cc(ptr:end,:);
    AA = AA(ptr:end,:);

    pars.fid = 0;
    [xx3,yy3,info3] = sedumi(AA',[],cc,KK,pars);
    feas3 = info3.pinf==0 && info3.dinf==0 && ...
        (info3.numerr==0 || info3.numerr==1);

    if feas1 && feas2 && feas3
        x0(ldv+1:end) = yy3(:);
    else
        warning('Initial point not feasible');
        x0  = [];
    end        
else
    x0 = [];
end

% ----- Solve BMI with PENBMI
[fopt,xopt,u,iflag,niter,feas] = bmi2pen(A0,A,K,fobj,msizes,x0);
bopt = xopt(1);
gopt = xopt(2);
s1opt = subs(s1,c1, xopt(2+(1:lc1))' );
s2opt = subs(s2,c2, xopt(2+lc1+(1:lc2))' );
Vopt = subs(V,cV, xopt(2+lc1+lc2+(1:lcV))' );

