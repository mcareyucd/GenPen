function [ beta_opt, C_opt, y_hat, var_beta, pars] = Gen_Pen_Non(argvals, y, basisobj, Term)
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
%pars         - Vector containing the optimal complexity parameter, Sum
%squared errors and degrees of freedom

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

%Non-Homogenous Term
K             = size(Phi,2);
V             = zeros(K,(nderiv+1));
Theta         = zeros(K.*(nderiv),1);
 for i = 1:(nderiv+1)
    V(:,i) = 1.75.*full(eval_basis(0.1, basisobj,i-1)); 
 end
Theta = reshape(V(:,1:nderiv),nderiv*K,1);

%Compute A and Tau matrices required for optimal parameters of the ODE
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
beta   = cell(100);
lambda = cell(100);
C      = cell(100);

%Set up the optimisation procedure
C{1} =  stinv((B'*B),B'*y);
Im   =  eye(nderiv);

%Compute the parameters for iteration
i   = 2;
eps = repmat(1e-6,nderiv,1);
con = 1;

while ( con && i<100 )
    
    %Compute beta
    Icrossc = kron(Im,C{i-1});
    beta{i} = Term.*(stinv(Icrossc'*(A)*Icrossc,Icrossc'*(-Tau*C{i-1}+Theta)));
    %Compute the penalty given beta
    R       = pen(beta{i},Psi,Ome,Phi);
    %Compute the optimal lambda using GCV
    lambda{i} = min_gcv_non(R,V,beta{i},B,y,Psi,Ome,A,Tau,Theta);
    %Compute the parameters of the smooth
    M       = B'*B + lambda{i}.*R;
    C{i}    = stinv((M),(B'*y+lambda{i}.*(V(:,1:nderiv)*beta{i}+V(:,nderiv))));
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


%% Partial derivatives for Variance estimates

%Partial derivative of c w.r.t. beta 
[c_beta] = dc_dbetaNH(C_opt,beta_opt,lambda_opt,Psi,Ome,M,V);
[beta_y] = dbeta_dyNH(beta_opt,A,Tau,y,stinv(M,B'),Theta);

S          = B*(c_beta*beta_y+stinv(M,B'));
df         = round(100*trace(eye(n)-S))./100;
pars       = [lambda_opt, SumSE, n-df]';

%% Variance Estimates

SigmaSq      = (SumSE./(n-df));
 
%Variance of the parameters of the ODE eqn (18) Parameter Estimation for differential equations paraper
%(2007) Ramsay, Hooker, Campbell, Cao.
var_beta = SigmaSq.*diag(stinv((B*c_beta)'*(B*c_beta),eye(nderiv))); 


end


