function [x]=stinv(A,b)
% Caluate the soltion to linesr system Ax=b, x=A-1b
%[x]=stinv(A,b)
%-----------Inputs--------------------------------
% A matrix to be inverted
% b vector
%-----------Outputs-------------------------------
% x vector x=A^(-1)b 

 if max(max(abs(A-A')))/max(max(abs(A))) > 1e-10
     error('First argument is not symmetric.');
 else
     A = (A + A')./2;
 end
%  Choleski decomposition Asym = Rmat'*Rmat;
[R,p] = chol(A); 
if p>0
    x=pinv(A)*b;    
else
    temp = R'\b;
    x = R\temp;
end

end
 

