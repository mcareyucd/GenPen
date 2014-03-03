function [beta_y]=dbeta_dy(beta_opt,A,Tau,y,MB)
%The partial derivative of the parameters of the ODE w.r.t y
%[beta_y]=dbeta_dy(beta_opt,A,Tau,y,MB)
%-----------Inputs--------------------------------
%beta_opt     - the parameters of the ODE
%S            - map of y to the coeffients of the smooth.
%A            - A matrix required for optimal parameters of the ODE
%Tau          - A matrix required for optimal parameters of the ODE
%y            - observaitions
%-----------Outputs-------------------------------
%beta_y       - The partial derivative of beta w.r.t y


%Set up preliminary matrices
 nderiv      = length(beta_opt);
 n           = length(y);
 beta_y      = zeros(nderiv,n);
 Im          = eye(nderiv);
% I kronceker c
 Icrossc     = kron(Im,MB*y);
% the inverted matrix in beta calulation
 cAc         = Icrossc'*A*Icrossc;

for i = 1:n
    %vector of zeros 
    e_n         = zeros(n,1);
    e_n(i)      = 1;
    % partial of I kroneker c w.r.t y
    Icrosse_n   = kron(Im,MB*e_n);
    % derivative of whats inside the inverse
    DcAc        = Icrosse_n'*A*Icrossc+Icrossc'*A*Icrosse_n;
    beta_y(:,i) = stinv(cAc,(-DcAc*beta_opt)+Icrosse_n'*(-Tau*MB*y)+Icrossc'*(-Tau*MB*e_n));
end
