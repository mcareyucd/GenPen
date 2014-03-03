function [optlambda]=min_gcv(R,B,y,beta_opt,Psi,Ome,Tau,A)% GCV as in Cao and Ramsay (2007)
%Computes the optimal lambda based on minimising generalised cross validation 
%[optlambda]=min_gcv(R,B,y)
%-----------Inputs--------------------------------
%R         - penalty matrix (Inner product of the ODE).
%B         - evaluated B-spline basis function.
%y         - data.
%-----------Outputs-------------------------------
%optlambda - optimal lambda

%Set of possible lambdas to select from
lambdas = 10.^(-6:1:6);
nlam    = length(lambdas);

%Calulate the geralised cross validation
gcv     = zeros(nlam,1);
Err     = zeros(nlam,1);
dfe     = zeros(nlam,1);
n       = length(y);

for j = 1:nlam
    %Evaluate M which maps y to c
    M       =  B'*B + lambdas(j).*R;
    MB      =  stinv(M,B');
    c       =  MB*y;
    %Evaluate S which maps y^hat to y
    [c_beta] = dc_dbeta(c,beta_opt,lambdas(j),Psi,Ome,M);
    [beta_y] = dbeta_dy(beta_opt,A,Tau,y,MB);
    S        = B*(c_beta*beta_y+MB); 
    %Compute the 1/n * Sum of squared Errors
    Err(j)    = (y-B*c)'*(y-B*c);
    %Compute the 1/n * error degrees of freedom
    dfe(j)     = n-round(100*trace(S))./100;
    %Compute the generalised cross validation
    gcv(j) = (Err(j)./(dfe(j).^2));
end

%Select the lambda that minimises gcv
optlambda = lambdas(gcv==min(gcv));
optlambda = optlambda(1);

end