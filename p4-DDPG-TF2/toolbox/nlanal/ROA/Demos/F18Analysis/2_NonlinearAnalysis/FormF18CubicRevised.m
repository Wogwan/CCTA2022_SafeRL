%
% This file creates a cubic approximation to the closed-loop F/A-18 full 
% nonlinear REVISED model. This is done in the following way: 
%
%   1.  Sample the f18clpfull.m within the operational limit of the states.
%   2.  This sampling evaluates the xcldot and return a vector evaluated
%   at the sampled value. 
%   3.  Subtract the linear terms and preserve linearity approximate the
%   higher order terms than linear. 
%   4.  Estimate  upto cbic approximation. Given xcldot @ x_c = x_0 xcldot
%   can be written as :
%     xcldot|x_0 = monomials(x_c)|x_0 * {coeff}'; coeff needs to be detemined
%    This is form as a LS problem for each state-derivative  
%    Forming the monomials(x_c) can be done in two ways: 
%         (a) Include all monomials
%         (b) Study the structure of the closed-loop and use the knowledge.
%         
% Abhijit 03/04/2010


% To ensure no variables from Form..Baseline gets carried over
clear all

% ------------- Cases for Trimming Nonlinear Model 
 
Level       = 0; 
CoordTurn   = 1;
UnCoordTurn = 0;
opts.poly   = 0; 

% ---------- Generate Trim Values / Open-loop linear model
GenF18LinModel; 

%----------- Generate Linear Closed-Loop Models / Load SS for Controller
GenF18LinCLPModel; 


% ---------- Revised Controller  : Build up Options 
opts.xtrim = xtrim;
opts.utrim = utrim; 
opts.Acl   = CLR.A; 
opts.K = CssR; 

% ---------- Form States in PVAR Format
pvar xcR
x = [beta ; alpha; p ;q ; r ;phi; xcR];

% ---------- Output matrix:
Jy = jacobian(y,u);
Hy = y-Jy*u;

% ---------- generate initial condition 
% XXX ----- Formulate bunch of ICs? Think about what is the appropiate way
% of sampling for this IC

Nvar = 7; 
Xdata.range = [-10 10; -5 45; -35 35; -30 30; -15 15; 10 60; -20 20]*d2r; 
Xdata.NPts =  [5; 6; 5; 5; 5; 4; 4];
% Xdata.NPts =  [5; 5; 5; 4; 5; 3; 3];
% Xdata.NPts =  [3; 3; 3; 3; 3; 3; 3];

X = cell(1,Nvar); 
for i1 = 1:Nvar
    X{:,i1} = linspace(Xdata.range(i1,1),Xdata.range(i1,2),Xdata.NPts(i1)); 
end

% --- Generate Gridded Data
Xgrid = cell(1,Nvar); 
[Xgrid{:}]=ndgrid(X{:});
Xdata.X = reshape(cat(length(Xgrid)+1,Xgrid{:}),[],length(Xgrid)) ;
Xdata.X = Xdata.X'; 

% -- total npts for grid
Xdata.Npts = prod(Xdata.NPts); 
        
xcldot2 = zeros(7,Xdata.Npts); 

for i1 = 1:(Xdata.Npts)
    
    x0 = Xdata.X(:,i1); 
    
    % ---------- Evaluate  output matrix and inverse of M 
    Jeval = double(subs(Jy,x,x0));
    Heval = double(subs(Hy,x,x0)); 
    Meval  = double(eye(7) + Jeval*opts.K.D); 
    Minv  = inv(Meval); 

    % ---------- Build up Options 

    opts.J = Jeval; 
    opts.H = Heval; 
    opts.Minv = Minv; 

    % ---------- Evaluate the Full Nonlinear Closed Loop Model 
    [xcldot2(:,i1),xcldot,f1,g1]  = f18Clpfull([],x0,opts);


end
%-----------------------------------------------------------------------
% Approximation Step 

% ----------- COMMENTS:
% ----------- Approximate the higher order terms xcldot2 with upto cubic
% degree polynomial .The Structure of nonlinearities has been investigated
% in the file F18NL_Structure.m file.  Run the file to understand how the
% what functional form of nonlinearities can be expected in each state
% deriavtive. See Notebook # 3 page 80

% ---------- Create Monomials to be used for the approximation 
a0   = monomials([alpha;beta;phi],2:3); 
a1   = monomials([alpha],2:3); 
a11  = monomials([alpha],1:2);
a2   = monomials([alpha;beta],2:3);
a22  = monomials([alpha;beta],1:2);
a33  =  monomials([phi],1:2); 

% % ---------- Create the range of data where LS to be fitted in 
% n = 15; 
betadata  = Xdata.X(1,:); 
alphadata = Xdata.X(2,:); 
pdata     = Xdata.X(3,:); 
qdata     = Xdata.X(4,:); 
rdata     = Xdata.X(5,:); 
phidata   = Xdata.X(6,:); 
xcRdata   = Xdata.X(7,:);


% ---------- Betadot 
% 
H1 = [a0;a22*p ;a22*xcR ;a22*q ;a22*r];
pbetadot = polydecvar('a',H1,'vec'); 
xbetadot = [beta;       alpha;      p;     q;      r;      phi;        xcR];
Xdatabeta = [betadata; alphadata; pdata;  qdata;  rdata;  phidata;  xcRdata];
Ydatabeta = xcldot2(1,:); 

[pbetafit,cbetafit,infobeta] = pdatafit(pbetadot,xbetadot,Xdatabeta,Ydatabeta);


% ---------- Alphadot 
% 
H2 = [a0;a22*p;a22*r;a22*q]; 
palphadot = polydecvar('b',H2,'vec'); 
xalphadot =  [beta;      alpha;     p;     q;      r;      phi]; 
Xdataalpha = [betadata;  alphadata; pdata; qdata; rdata;  phidata];
Ydataalpha = xcldot2(2,:); 

[palphafit,calphafit,infoalpha] = pdatafit(palphadot,xalphadot,Xdataalpha,Ydataalpha);

% ---------- pdot 
% 
H3 = [a2; p*q ; q*r; a11*p; a11*r; a11*xcR ; a11*phi];
ppdot = polydecvar('c', H3,'vec'); 
xpdot =  [beta;      alpha;     p;      q;      r;      phi;     xcR]; 
Xdatap = [betadata; alphadata; pdata;  qdata;   rdata;  phidata; xcRdata];
Ydatap = xcldot2(3,:);

[ppfit,cpfit,infop] = pdatafit(ppdot,xpdot,Xdatap,Ydatap);


% ---------- qdot 
% 
H4 = [a1; a11*q;  p*r;  r^2;  p^2]; 
pqdot = polydecvar('d',H4,'vec'); 
xqdot = [alpha;p;q;r]; 
Xdataq = [alphadata;pdata;qdata;rdata];
Ydataq = xcldot2(4,:);

[pqfit,cqfit,infoq] = pdatafit(pqdot,xqdot,Xdataq,Ydataq);

% ---------- rdot 
% 
H5 = [a2 ; p*q ; q*r; a11*p; a11*r; a11*xcR]; 
prdot = polydecvar('c', H5,'vec'); 
xrdot =  [beta;      alpha;     p;      q;      r;      phi;     xcR]; 
Xdatar = [betadata; alphadata; pdata;  qdata;   rdata;  phidata; xcRdata];
Ydatar = xcldot2(5,:);

[prfit,crfit,infor] = pdatafit(prdot,xrdot,Xdatar,Ydatar);



% ---------- phidot 
% 
H6 = [a33*r; a33*q];
pphidot = polydecvar('e', H6,'vec'); 
xphidot = [phi;q;r]; 
Xdataphi = [phidata;qdata;rdata];
Ydataphi = xcldot2(6,:);

[pphifit,cphifit,infophi] = pdatafit(pphidot,xphidot,Xdataphi,Ydataphi);


% -------------------------------------------------------------------------
% Form the revised model 

xcldotR = [pbetafit; palphafit; ppfit; pqfit ;  prfit; pphifit; 0] + CLR.A*x;


% --------------- Save Data
V2 = strcat('2_Mar14_PolyF18RevisedModel_Phi_',num2str(phi0*r2d));
save(V2,'xcldotR','CLR','xtrim','utrim','Pss','Xdata','x')