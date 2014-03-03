function [ beta_opt, C_opt, y_hat, var_beta, var_c, pars] = Gen_Pen(argvals, y, basisobj, Term)
%Computes the optimal beta from joint optimisation routine described in ODE paper.
%[beta] = generalised_smooth(argvals, y, fdParobj)
%-----------Inputs--------------------------------
%argvals      - the locations.
%y            - the observations.
%fdParobj     - functional data parameter object.
%Term         - Boolean vector indicating the derivatives in ODE model.
%-----------Outputs-------------------------------
%beta_opt     - Optimal parameters of the ODE
%C_opt        - Optimal paramters of the smooth
%y_hat        - The estimated functional entity
%var_beta     - Variance of the Optimal parameters of the ODE
%var_c        - Variance of the Optimal parameters of the smooth
%var_f        - Variance of the functional entity
%pars         - Vector containing the optimal complexity parameter

% if nargin < 4
% error('There is not at least four arguments.');
% end

%  check ARGVALS
[argvals, n] = argcheck(argvals);
rng(1)       = min(argvals);
rng(2)       = max(argvals);

%  check Y
[y]          = ycheck(y, n);

%  extract the ORDER of the ODE
nderiv      = length(Term);

%  evaluate the B-SPLINE
B           = eval_basis(argvals, basisobj);
B           = full(B);

%  compute products of derivatives 0 to nderiv-1 (Psi), 
%          products with D^nderiv (Ome)
%          and inner product of D^nderiv (Phi)
[Psi,Ome,Phi] = bsplinepenJ(basisobj, nderiv, rng);

%Compute A and Tau matrices required for optimal parameters of the ODE
K             = size(Phi,2);
Ik            = eye(K);
A             = zeros(K.*(nderiv),K.*(nderiv));
Tau           = zeros(K.*(nderiv),K);

for i = 1:nderiv
    if(Term(i)==1)
        e_i       = zeros(nderiv,1);
        e_i(i)    = 1;
        Icrosse_i = kron(Ik,e_i);
        for j=1:nderiv
            if(Term(j)==1)
                e_j        = zeros(nderiv,1);
                e_j(j)     = 1;
                Icrosse_j  = kron(Ik,e_j);
                A(((i-1)*(K)+1):((i)*(K)),((j-1)*(K)+1):((j)*(K))) = Icrosse_i'*Psi*Icrosse_j + Icrosse_j'*Psi*Icrosse_i;
            end
        end
        Tau(((i-1)*(K)+1):((i)*(K)),:) = Ome'*Icrosse_i;
    end
end

%% Joint Optimiation 

%Compute the parameters for iteration
i      = 2;
eps    = repmat(1e-5,nderiv,1);
con    = 1;
beta   = cell(100);
lambda = cell(100);
C      = cell(100);

%Set up the optimisation procedure
Im   =  eye(nderiv);
C{1} =  stinv((B'*B),B'*y);

while ( con && i<100 )
    
    %Compute beta
    Icrossc = kron(Im,C{i-1});
    beta{i} = Term.*(stinv(Icrossc'*(A)*Icrossc,-Icrossc'*Tau*C{i-1}));
    %Compute the penalty given beta
    R       = pen(beta{i},Psi,Ome,Phi);
    %Compute the optimal lambda using GCV
    %lambda{i} = GCV(R,B,y,beta{i},Psi,Ome,Tau,A,K,n,nderiv);
    lambda{i} = min_gcv(R,B,y,beta{i},Psi,Ome,Tau,A);
    %Compute the parameters of the smooth
    M       = B'*B + lambda{i}.*R;
    C{i}    = stinv((M),(B'*y));
    %Evaluate the accuarcy of the estimates 
    if(i>2)
    dif     = abs(beta{i}-beta{i-1});
    con     = sum(dif>eps)~=0;
    end
    i       = i+1;

end

%Save optimal parameters
opt        = i-1;
beta_opt   = beta{opt};
C_opt      = C{opt};
lambda_opt = lambda{opt};
y_hat      = B*C_opt;
SumSE      = (y-y_hat)'*(y-y_hat);
R          = pen(beta_opt,Psi,Ome,Phi);
M          = (B'*B + (lambda_opt).*R);
MB         = stinv(M,B');

%% Partial derivatives for Variance estimates
ind      = find(Term==1);

%Partial derivative of c w.r.t. beta 
[c_beta] = dc_dbeta(C_opt,beta_opt,lambda_opt,Psi,Ome,M);
c_beta   = c_beta(:,ind);

%Partial derivative of beta w.r.t. y 
[beta_y] = dbeta_dy(beta_opt,A,Tau,y,MB);
beta_y   = beta_y(ind,:);

S          = B*(c_beta*beta_y+MB);
df         = trace(S);
dfe        = n-round(100*trace(S))./100;
pars       = [lambda_opt, SumSE, df]';


%% Variance Estimates
SigmaSq      = (SumSE./(dfe));

%Variance of the parameters of the ODE eqn (18) Parameter Estimation for differential equations paraper
%(2007) Ramsay, Hooker, Campbell, Cao.
var_beta = SigmaSq.*diag(stinv((B*c_beta)'*(B*c_beta),eye(nderiv))); 

var_c = SigmaSq.*diag((c_beta*beta_y+MB)*(c_beta*beta_y+MB)');

end


