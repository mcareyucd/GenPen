function [Ps,Ome,Phi] = bsplinepen(basisobj, nderiv, rng)

%BSPLINEPEN computes the bspline penalty matrix for penalty LFDOBJ.
%  Arguments:
%  BASISOBJ  ... a basis object
%  LFDOBJ    ... a linear differential operator object.  
%                Default int2Lfd(2)
%  RNG       ... A range over which the product is evaluated
%  SPARSEWRD ... if 1, return penaltymatrix in sparse storage mode.

%  Last modified:  1 November 2007

%  check BASISOBJ

if ~isa_basis(basisobj)
    error('BASISOBJ is not a basis object.');
end

%  check basis type

type = getbasistype(basisobj);
if ~strcmp(type, 'bspline')
          error('basisobj not of type bspline');
end

%  set up default value for RNG

range = getbasisrange(basisobj);
if nargin < 3, rng = range;  end

%  get basis information

nbasis  = getnbasis(basisobj);
params  = getbasispar(basisobj);

%  if there are no internal knots, use the monomial penalty

if isempty(params)
          basisobj   = create_monomial_basis(range, nbasis, 0:nbasis-1);
          penaltymat = monompen(basisobj, Lfdobj, rng);
          return;
end

%  normal case:  PARAMS is not empty

breaks    = [rng(1),params,rng(2)];  %  break points
nbreaks   = length(breaks);
norder = nbasis - length(params);

%  check break values

if length(breaks) < 2
          error('The length of argument breaks is less than 2.');
end

%  Set up the knot sequence
knots =[range(1)*ones(1,norder),params,range(2)*ones(1,norder)];

%Set up blank matrices for the outputs
Psi = zeros(nbasis*(nderiv),nbasis*(nderiv));
Ome_a = zeros(nbasis,nbasis*(nderiv));
Phi = zeros(nbasis,nbasis);
an = zeros((nderiv+1)*(nderiv+1),1);

for k = 1:nbasis
    %coefficients of the piece-wise polynomials
    %B(k,t)(x) =
    %   COEFF(i,1)*(x-T(i))^{k-1} + ... + COEFF(i,k-1)*(x-T(i)) +
    %   COEFF(i,k), for x in [T(knots[i]),T(knots[i+1])).
    Coeff1 = ppBspline(knots(k:(k+norder)));
    for j = k:(k+norder-1)
        if(j<=nbasis)
            %coefficients of the piece-wise polynomials
            Coeff2 = ppBspline(knots(j:(j+norder)));
            % evaluate all products of the derivatives
            an = sintpoly3(Coeff1,Coeff2,knots,k,j,nderiv);
            Psi_ans = zeros((nderiv),(nderiv));
            Ome_ans = zeros(1,(nderiv));
            for i=0:(nderiv-1)
                Psi_ans(i+1,:)=an((i*(nderiv+1)+1):(i*(nderiv+1)+nderiv));
                Ome_ans(:,i+1)=an((i*(nderiv+1)+(nderiv+1)))+an(((nderiv+1)*(nderiv+1)-(nderiv+1)+i+1));
            end
            Psi( (k-1)*nderiv+(0:(nderiv-1))+1, (j-1)*nderiv+(0:(nderiv-1))+1 ) = Psi_ans;
            Ome_a(k, (j-1)*nderiv+(0:(nderiv-1))+1 ) = Ome_ans;
            Phi(k,j)=an((nderiv+1)*(nderiv+1));
        end
    end
end

 %make the matrices symetric
  Phi=Phi + triu(Phi, 1)';
  
  Ome = zeros(nbasis*(nderiv),nbasis);
  Ps  = zeros(nbasis*(nderiv),nbasis*(nderiv));
  
  for i=1:nderiv
      Om_1=Ome_a(:,i:nderiv:((nbasis)*nderiv));
      Om_1=Om_1 + triu(Om_1,1)';
      Ome(i:nderiv:((nbasis)*nderiv),:)=Om_1;
      for j=1:nderiv
          Ps_1=Psi(i:2:(nbasis*nderiv),j:2:(nbasis*nderiv));
          Ps_1=Ps_1 + (-1)^(i+j)*triu(Ps_1,1)';
          Ps(i:2:(nbasis*nderiv),j:2:(nbasis*nderiv))=Ps_1;
      end
  end

  end



