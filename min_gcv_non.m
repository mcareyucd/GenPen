function [optlambda]=min_gcv_non(R,V,beta,B,y,Psi,Ome,A,Tau,Theta)
% GCV as in Cao and Ramsay (2007)
%Computes the optimal lambda based on minimising generalised cross validation 
%[optlambda]=min_gcv(R,B,y)
%-----------Inputs--------------------------------
%R         - penalty matrix (Inner product of the ODE).
%B         - evaluated B-spline basis function.
%y         - data.
%-----------Outputs-------------------------------
%optlambda - optimal lambda

%Set of possible lambdas to select from
lambdas = 10.^(-3.5:1:-1);
nlam    = length(lambdas);    

%Calulate the geralised cross validation
gcv     = zeros(nlam,1);
Err     = zeros(nlam,1);
dfe     = zeros(nlam,1);
n       = length(y);
nderiv  = length(beta);

for j = 1:nlam
    %Evaluate M which maps y to c
    M      = B'*B + lambdas(j).*R;
    C      = stinv(M,B'*y+lambdas(j).*(V(:,1:nderiv)*beta+V(:,nderiv)));
    %Evaluate S which maps y^hat to y
    F      = B*C; 
    %Compute the 1/n * Sum of squared Errors
    Err(j)    = (y-F)'*(y-F);
    %Compute the 1/n * error degrees of freedom
    %Partial derivative of c w.r.t. beta 
    [c_beta] = dc_dbetaNH(C,beta,lambdas(j),Psi,Ome,M,V);
    [beta_y] = dbeta_dyNH(beta,A,Tau,y,stinv(M,(B')),Theta);
    S      = B*(c_beta*beta_y+stinv(M,B')); 
    dfe(j) = round(100*trace(eye(n)-S))./100;
    %Compute the generalised cross validation
    gcv(j) = n.*(Err(j)./(dfe(j).^2));
end
% figure();
% plot(log(lambdas),log(gcv),'b')
% hold on;
% plot(log(lambdas),log(Err),'-r');
% plot(log(lambdas),log(dfe),'-g');
% pause;

%Select the lambda that minimises gcv
optlambda = lambdas(gcv==min(gcv));
optlambda = optlambda(1);

end