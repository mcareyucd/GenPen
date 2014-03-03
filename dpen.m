function [dR]=dpen(beta,e_j,Psi,Ome)
%evaluates the derivative of the penalty w.r.t beta
%[R]=pen(beta,Psi,Ome,Phi)
%beta      - current value of beta.
%Psi       - products of derivatives 0 to nderiv-1.
%Ome       - products of derivatives nderiv with derivatives from 0 to nderiv-1.
%Phi       - inner product of nderiv.
%B         - evaluated B-spline basis function.
%y         - data.

%  compute roughness penalty matrix
    Id              = eye(size(Ome,2));
    Icrossbeta      = kron(Id,beta);
    Icrossej        = kron(Id,e_j);
    Psi_Icrossbeta  = Psi*Icrossbeta;
    W               = 2.*Psi_Icrossbeta + Ome;
    dR              = Icrossej'*W;    
    
 end