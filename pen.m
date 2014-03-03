function [R]=pen(beta,Psi,Ome,Phi)
%evaluates the penalty
%[R]=pen(beta,Psi,Ome,Phi)
%beta      - current value of beta.
%Psi       - products of derivatives 0 to nderiv-1.
%Ome       - products of derivatives nderiv with derivatives from 0 to nderiv-1.
%Phi       - inner product of nderiv.
%B         - evaluated B-spline basis function.
%y         - data.

%  compute roughness penalty matrix
    Id              = eye(size(Phi,2));
    Icrossbeta      = kron(Id,beta);
    Psi_Icrossbeta  = Psi*Icrossbeta;
    W               = Psi_Icrossbeta + Ome;
    R               = Icrossbeta'*W  + Phi;    
    
 end