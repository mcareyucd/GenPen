function [c_beta]=dc_dbeta(C_opt,beta_opt,lambda_opt,Psi,Ome,M)
%The partial derivative of the parameters of smooth w.r.t the parameters of
%the ODE
%[c_beta]=dc_dbeta(C_opt,beta_opt,lambda_opt,Psi,Ome,M)
%-----------Inputs--------------------------------
%C_opt        - the parameters of the smooth
%beta_opt     - the parameters of the ODE
%lambda_opt   - the complexity parameter
%Psi          - products of derivatives 0 to nderiv-1
%Ome          - products of derivatives 0 to nderiv-1 with D^nderiv
%M            - A matrix required for optimal parameters of the ODE
%y            - observaitions
%-----------Outputs-------------------------------
%c_beta       - The partial derivative of c w.r.t beta

nderiv  = size(beta_opt,1);
K       = size(C_opt,1);
c_beta  = zeros(K,nderiv);

for j = 1:nderiv
    e_j            = zeros(nderiv,1);
    e_j(j)         = 1;
    dR_dbeta       = dpen(beta_opt,e_j,Psi,Ome);
    c_beta(:,j)    = -lambda_opt.*stinv(M,(dR_dbeta'*C_opt));
end