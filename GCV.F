#include "fintrf.h"

C======================================================================

C======================================================================
C     Functions
      subroutine stinv(A, B, In, nA, nB)

      !External
      integer*8 nA, nb
      real*8 A(1:nA,1:nA), B(1:nA,1:nb)
      real*8 In(1:nA,1:nb)
      

      !Internals
      integer*8 i, j, INFO, lmax
      integer*8 one, rank, lwork, IPIV(nA)
      real*8 work(1000), Ain(1:nA,1:nA)
      real*8 S(1:2*nA), rcond, Bin(1:nA,1:nb) 
      real*8 CPA(1:nA,1:nA), WR(1:nA), WI(1:nA), D(1,1)
      real*8 R(1:nA,1:nA), Temp(1:nA,1:nb), tR(1:nA,1:nA)
      logical p
      
      In = 0
      
      CPA = A
      Ain = A
      Bin = B
      one = 1
      lwork = -1
      
      call dgeev('N','N',nA,CPA,nA,WR,WI,D,one,D,one,work,lwork,INFO)
      
      lmax = 1000
      lwork = min(lmax, int( work( 1 ) ) )
      
      call dgeev('N','N',nA,CPA,nA,WR,WI,D,one,D,one,work,lwork,INFO)
     
      p = .TRUE.
      do i=1,nA 
      if (WR(i) .LT. 0.0) then
      p = .FALSE.
      end if
      end do
      
      In = 0
      
      if( p ) then
      call DPOTRF( 'U', nA, A, nA, INFO )
      R=0
      do i=1,nA
      do j=i,nA
      R(i,j) = A(i,j)
      end do
      end do 
      tR = transpose(R)
      call DGESV( nA, nb, tR, nA, IPIV, B, nA, INFO )
      Temp = 0
      Temp = B
      call DGESV( nA, nb, R, nA, IPIV, Temp, nA, INFO )
      In = Temp
      else
      rank  = 0.01d0
      call dgelss(nA,nA,nB,A,nA,B,nA,S,rcond,rank,work,lwork,INFO)
      In = B
      end if 
      
      A=Ain
      B=Bin
      
      return

      end subroutine stinv

      subroutine kron(KAb, A, b, nA, nb)

      integer*8 i, j, nA, nb
      real*8 A(1:nA,1:nA), KAb(1:(nb*nA),1:nA) 
      real*8 b(1:nb)

      do i=1,nA
      do j=1,nA
        KAb((nb*(i-1)+1):(nb*(i-1))+nb,j) = A(i,j)*b(1:nb)
      end do
      end do
      
      return

      end subroutine kron

      subroutine iden(n,A)

      integer*8 n, i, j
      real*8 A(1:n,1:n)

      do i =1,n
      do j =1,n
      A(i,j)=0
      if(i.eq.j) then
      A(i,j)=1
      end if
      end do
      end do

      return
      end subroutine iden

      subroutine trace(n,A,b)

      integer*8 n, i
      real*8 A(1:n,1:n), b
      b=0
      do i =1,n
      b=b+A(i,i)
      end do

      return
      end subroutine trace
    
      subroutine dpen(beta,e_j,Psi,Ome,dR,nderiv,K)
      !evaluates the derivative of the penalty w.r.t beta
      ![dR]=dpen(beta,Psi,Ome,Phi)
      !beta      - current value of beta.
      !Psi       - products of derivatives 0 to nderiv-1.
      !Ome       - products of derivatives nderiv with derivatives from 0 to nderiv-1.
      !Phi       - inner product of nderiv.
      !B         - evaluated B-spline basis function.
      !y         - data.

      !  compute roughness penalty matrix
      integer*8 nderiv, K

      real*8 beta(1:nderiv), e_j(1:nderiv)
      real*8 dR(1:K,1:K)
      real*8 Psi(1:K*nderiv,1:K*nderiv)
      real*8 Ome(1:K*nderiv,1:K)

      real*8 idene(1:K,1:K)
      real*8 Icrossbeta(1:K*nderiv,1:K)
      real*8 Icrossej(1:K*nderiv,1:K)
      real*8 W(1:K*nderiv,1:K)
      real*8 Psi_Icrossbeta(1:K*nderiv,1:K)

      call  iden(K,idene)
      call  kron(Icrossbeta,idene,beta,K,nderiv)
      call  kron(Icrossej,idene,e_j,K,nderiv)

      Psi_Icrossbeta  = matmul(Psi,Icrossbeta)
      W               = 2.d0*Psi_Icrossbeta + Ome;
      dR              = matmul(transpose(Icrossej),W)

      return

      end subroutine dpen

      subroutine dbeta_dy(beta_opt,A,Tau,y,MB,beta_y,K,n,nderiv)

      integer*8 K, n, nderiv, i, j, one

      real*8 Tau(1:K*nderiv,1:K) 
      real*8 A(1:K*nderiv,1:K*nderiv)
      real*8 y(1:n), beta_opt(1:nderiv)
      real*8 beta_y(1:nderiv,1:n)
      real*8 MB(1:K,1:n)

      real*8 Im(1:nderiv,1:nderiv), c(1:K)
      real*8 Icrossc(1:nderiv*K,1:nderiv)
      real*8 cAc(1:nderiv,1:nderiv)
      real*8 Icrosse_n(1:nderiv*K,1:nderiv)
      real*8 en(1:n), IcAc(1:nderiv,1:nderiv)
      real*8 IenA(1:nderiv,1:K*nderiv) 
      real*8 ST(1:nderiv), TT(1:nderiv)
      real*8 FT(1:nderiv), DcAc(1:nderiv,1:nderiv)
      real*8 cA(1:nderiv,1:K*nderiv)
      real*8 AIen(1:K*nderiv,1:nderiv)
      real*8 MBen(1:K), STer(1:nderiv,1)
      real*8 In(1:nderiv,1)

C      Set up preliminary matrices
      call  iden(nderiv,Im)
C     I kronceker c
      c = matmul(MB,y)
      call kron(Icrossc,Im,c,nderiv,K)
C     the inverted matrix in beta calulation
      cA  = matmul(transpose(Icrossc),A)
      cAc = matmul(cA,Icrossc)

      one = 1
      
      do i = 1,n
      en  = 0
      en(i) = 1.d0
      MBen  = matmul(MB,en)
C     partial of I kroneker c w.r.t y
      call kron(Icrosse_n,Im,MBen,nderiv,K)
C     derivative of whats insidene the inverse
      IenA  = matmul(transpose(Icrosse_n),A)
      AIen  = matmul(A,Icrosse_n)
      DcAc  = matmul(IenA,Icrossc)+matmul(transpose(Icrossc),AIen)

      ST = matmul(-DcAc,beta_opt)
      TT = matmul(transpose(Icrosse_n),matmul(-Tau,c))
      FT = matmul(transpose(Icrossc),matmul(-Tau,MBen))

      STer(1:nderiv,1) = ST(1:nderiv)+TT(1:nderiv)+FT(1:nderiv)
      
      call stinv(cAc, STer, In, nderiv, one)
     
      beta_y(:,i) = In(:,1)
     
      end do

      return

      end subroutine dbeta_dy

      subroutine dc_dbeta(pars,Psi,Ome,M,c_beta,K,n,nderiv)
      !The partial derivative of the parameters of smooth w.r.t the parameters of
      !the ODE
      ![c_beta]=dc_dbeta(C_opt,beta_opt,lambda_opt,Psi,Ome,M)
      !-----------Inputs--------------------------------
      !C_opt        - the parameters of the smooth
      !beta_opt     - the parameters of the ODE
      !lambda_opt   - the complexity parameter
      !Psi          - products of derivatives 0 to nderiv-1
      !Ome          - products of derivatives 0 to nderiv-1 with D^nderiv
      !M            - A matrix required for optimal parameters of the ODE
      !y            - observaitions
      !-----------Outputs-------------------------------
      !c_beta       - The partial derivative of c w.r.t beta

      integer*8 K, n, nderiv, j, one

      real*8 Psi(1:K*nderiv,1:K*nderiv)
      real*8 Ome(1:K*nderiv,1:K)
      real*8 M(1:K,1:K)
      real*8 c_beta(1:K,1:nderiv)
      real*8 y(1:n), beta_opt(1:nderiv) 
      real*8 C_opt(1:K), lambda_opt
      real*8 e_j(1:nderiv), In(1:K,1) 
      real*8 dR_dbeta(1:K,1:K), dRc(1:K,1)
      real*8 pars(1:K+nderiv+2)

      C_opt(1:K) = pars(1:K)
      beta_opt(1:nderiv) = pars((K+1):K+1+nderiv)
      lambda_opt = pars(K+nderiv+2)
      one=1
      
      do j = 1,nderiv
      e_j    = 0
      e_j(j) = 1
      call dpen(beta_opt,e_j,Psi,Ome,dR_dbeta,nderiv,K)
      dRc(1:K,1) = matmul(transpose(dR_dbeta),C_opt)
   
      call stinv(M, dRc, In, K, one)
      
      c_beta(:,j) = -lambda_opt*In(:,1)
      
      end do

      return

      end subroutine dc_dbeta

C     Computational subroutine

      subroutine GCV(R,B,y,beta_opt,Psi,Ome,Tau,A,optlambda,K,n,nderiv)

      integer*8 K, n, nderiv, i
      integer*4 nn, j, INFO

      real*8 R(1:K,1:K), B(1:n,1:K)
      real*8 Psi(1:K*nderiv,1:K*nderiv)
      real*8 Ome(1:K*nderiv,1:K)
      real*8 Tau(1:K*nderiv,1:K)
      real*8 A(1:K*nderiv,1:K*nderiv)
      real*8 y(1:n), beta_opt(1:nderiv)
      real*8 optlambda(1,1)

C     Interals
      real*8 lambdas(1:13), gcvv(1:13), c(1:K)
      real*8 c_beta(1:K,1:nderiv), res(1:n)
      real*8 M(1:K,1:K), IM(1:K,1:K) 
      real*8 MB(1:K,1:n), S(1:n,1:n)
      real*8 beta_y(1:nderiv,1:n), trS
      real*8 pars(1:K+nderiv+2), Tb(1:K,1:n)
      real*8 In(1:K,1:K)
      real*8 er(1:13), dfe(1:13)
      
      nn=n
      pars(K+1:K+1+nderiv) = beta_opt(1:nderiv)
      Tb =transpose(B)
      gcvv = 0
      dfe = 0
      er = 0
      
      do  j = 1,13 
       
      lambdas(j) = 10.d0**(j-7)
C     Evaluate M which maps y to c
      M       =  matmul(Tb,B) + lambdas(j)*R 
      call stinv(M, Tb, MB, K, n)
      c       = matmul(MB,y)
      
C     Evaluate S which maps y^hat to y
      pars(1:K)  = c
      pars(K+nderiv+2) = lambdas(j)
      call     dc_dbeta(pars,Psi,Ome,M,c_beta,K,n,nderiv)
      call     dbeta_dy(beta_opt,A,Tau,y,MB,beta_y,K,n,nderiv)
      S        = matmul(B,matmul(c_beta,beta_y)+MB)
      
C     Compute the 1/n * Sum of squared Errors
      res   = y-matmul(B,c)
      er(j) = dot_product(res,res)
      
C     Compute the 1/n * error degrees of freedom
      call    trace(n,S,trS)
      dfe(j)  = n-(anint(100.d0*trS)/100.d0)
      
C     Compute the generalised cross validenation
      gcvv(j) = (er(j)/(dfe(1)**2.d0))
      
      end do
      
      optlambda(1,1) = lambdas(minloc(gcvv,1))
      !optlambda = beta_y
      
      return

      end subroutine gcv


C-----------------------------------------------------------------------
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      mwPointer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetM, mxGetN
      integer   mxIsNumeric

C     Pointers to input/output mxArrays:
      mwPointer R_pr, B_pr, y_pr, beta_opt_pr, Psi_pr, Ome_pr, Tau_pr
      mwPointer A_pr, optlambda_pr

C     Array information:
      integer*8 K, n, nderiv
     
      !real*8 R(K*K), B(n*K), Psi(K*nderiv,K*nderiv), Ome(K*nderiv,K)
      !real*8 Tau(K*nderiv,K), A(K*nderiv,K*nderiv)
      !real*8 y(n), beta_opt(nderiv), optlambd(1)

      real*8 R(300,300), B(300,300), Psi(300,300)
      real*8 Ome(300,300), Tau(300,300), A(300,300)
      real*8 y(300,1), beta_opt(300,1), optlambda(1,1)
      !real*8,dimension(:,:),allocatable :: optlambda

C-----------------------------------------------------------------------
C     Check for proper number of arguments.
C      if (nrhs .ne. 3) then
C         call mexErrMsgidenAndTxt ('MATLAB:convhp:nInput',
C     +                           'Three inputs required.')
C      elseif (nlhs .gt. 1) then
C         call mexErrMsgidenAndTxt ('MATLAB:convhp:nOutput',
C     +                           'One output required.')
C      endif

C     Validenate inputs
      K       = mxGetM(prhs(1))
      n       = mxGetM(prhs(2))
      nderiv  = mxGetM(prhs(4))
      
C     Check that inputs are row vectors.
C      if(mx .ne. 1 .or. my .ne. 1) then
C         call mexErrMsgidenAndTxt ('MATLAB:convec:NonRowVector',
C     +                           'Inputs must be row vectors.')
C     Check size of the inputs.
C      elseif(nx .gt. 100 .or. ny .gt. 100) then
C         call mexErrMsgidenAndTxt ('MATLAB:convec:MaxElements',
C     +                 'Inputs must have less than 100 elements.')
C     Check if inputs are numeric.
C      elseif ((mxIsNumeric(prhs(1)) .ne. 1) .or.
C     +        (mxIsNumeric(prhs(2)) .ne. 1)) then
C         call mexErrMsgidenAndTxt ('MATLAB:convec:NonNumeric',
C     +                           'Inputs must be numeric.')
C      endif

C     Create the output array.
      plhs(1)     = mxCreateDoubleMatrix(1,1,0)
      R_pr        = mxGetPr(prhs(1))
      B_pr        = mxGetPr(prhs(2))
      y_pr        = mxGetPr(prhs(3))
      beta_opt_pr = mxGetPr(prhs(4))
      Psi_pr      = mxGetPr(prhs(5))
      Ome_pr      = mxGetPr(prhs(6))
      Tau_pr      = mxGetPr(prhs(7))
      A_pr        = mxGetPr(prhs(8))
      optlambda_pr = mxGetPr(plhs(1))

C     Load the data into Fortran arrays.
      call mxCopyPtrToReal8(R_pr,R,K*K)
      call mxCopyPtrToReal8(B_pr,B,n*K)
      call mxCopyPtrToReal8(y_pr,y,n)
      call mxCopyPtrToReal8(beta_opt_pr,beta_opt,nderiv)
      call mxCopyPtrToReal8(Psi_pr,Psi,K*nderiv*K*nderiv)
      call mxCopyPtrToReal8(Ome_pr,Ome,K*nderiv*K)
      call mxCopyPtrToReal8(Tau_pr,Tau,K*nderiv*K)
      call mxCopyPtrToReal8(A_pr,A,K*K*nderiv*nderiv)

C     Call the computational subroutine.
      call GCV(R,B,y,beta_opt,Psi,Ome,Tau,A,optlambda,K,n,nderiv)

C     Load the output into a MATLAB array.
      call mxCopyReal8ToPtr(optlambda,optlambda_pr,1)

      return
      end