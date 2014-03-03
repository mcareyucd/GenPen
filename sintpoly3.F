#include "fintrf.h"

C======================================================================
C
#if 0
C     
C     convh.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C
C     Convolves two input vectors
C      
C======================================================================
C     Functions
      subroutine fact(n,fac)

      integer*8 n,i,fac
      fac=1

      if(n>0) then

      do 10 i=1,n

      fac=fac*i

  10  continue

      endif

      return

      end

C     Computational subroutine

      subroutine sintpoly3(x,y,kn,k,j,nd,z,nx,ny,mx,my,nkn)

      integer*8 nx, ny, mx, my, nkn, ik, ij
      integer*8 t, t1, t2, i, r, s, p, q
      real*8  x(mx*nx), y(my*ny), kn(nkn), dkn(1)
      real*8  k(1), j(1), nd(1), z((int(nd(1))+1)*(int(nd(1))+1))
      integer*8 num, denom, der, nbasis

C     internal
      real*16 xx(mx*nx), yy(my*ny), kns(nx+ny-1)
      real*16 den(nx+ny-1), zx(nx+ny-1)
      real*16 drx(nx), dsy(ny)
      real*16 smx(1), ans((int(nd(1))+1)*(int(nd(1))+1))
      logical ind1, ind2

      xx=x
      yy=y

      ik=k(1)
      ij=j(1)
      der = nd(1)
      nbasis = nkn-nx

      do 10   q=1,((der+1)*(der+1))
         ans(q) = 0.00
         z(q)   = 0.00
  10     continue

      do 20 r=0,der

      do 30 s=0,der

      do 40 t=1,mx

      ind1 = (ik+t-ij).GT.0  .AND. ik.GE.nx
      ind2 = (ij-nx+t).LE.ik .AND. ik.LT.nx

      if ( ind1 .OR. ind2 ) then
     
      if ( ik.GE.nx .AND. ik.LE.nbasis) then

      dkn(1) = kn(ik+t)-kn(ik+t-1)
      t1 = t
      t2 = (ik+t-ij)

      else if (ik.LT.nx) then

      if (ij.LE.nx) then
      t1  = t
      t2  = t
      dkn(1) = kn(nx+t)-kn(nx+t-1)
      else
      t1  = ij-nx+t
      t2  = t
      dkn(1) = kn(nx+t1)-kn(nx+t1-1)
      endif

      else if (ik.GT.nbasis) then

      dkn(1) = kn(ik+t)-kn(nbasis)
      t1 = t
      t2 = (ik+t-ij)

      endif

      do 50 i=1,nx
        drx(i) = 0.00
        dsy(i) = 0.00
  50     continue

      do 60 i=1,(nx-r)
       call fact((nx-i),num)
       call fact((nx-i-r),denom)
       drx(r+i) = (num/denom)* xx((i-1)*mx+t1)
  60     continue

      do 70 i=1,(nx-s)
      call fact((nx-i),num)
      call fact((nx-i-s),denom)
      dsy(s+i) = (num/denom)* yy((i-1)*my+t2)
  70     continue

C     Initialize the output array
      do 80 i=1,(nx+ny-1)
         zx(i)  = 0
         den(i) = nx+ny-i
         kns(i) = (dkn(1))**(den(i))
  80    continue

C     Compute the convoluted vector
      do 100 i=1,nx
         do 90 p=1,ny
            zx(i+p-1) = zx(i+p-1) + drx(i) * dsy(p)
  90   continue
  100   continue

      smx(1) = 0.00
C     Integrate
      do 110 i=1,(nx+ny-1)
           smx(1) = smx(1)+(zx(i)*real(1,
     +       kind=selected_real_kind(4*precision(1.0)))/den(i))*kns(i)
  110   continue

           ans(r*(der+1)+s+1)=ans(r*(der+1)+s+1)*real(1,
     +       kind=selected_real_kind(4*precision(1.0)))+smx(1)*real(1,
     +       kind=selected_real_kind(4*precision(1.0)))

      endif

  40   continue

  30   continue

  20   continue

      z=ans

      return

      end

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
      mwPointer x_pr, y_pr, kn_pr, k_pr, j_pr, nd_pr, z_pr

C     Array information:
      integer*8 mx, nx, my, ny, nkn
      integer*8 one

      real*8 x(1000), y(1000), kn(1000), k(1), j(1), nd(1), z(100)

C-----------------------------------------------------------------------
C     Check for proper number of arguments.
C      if (nrhs .ne. 3) then
C         call mexErrMsgIdAndTxt ('MATLAB:convhp:nInput',
C     +                           'Three inputs required.')
C      elseif (nlhs .gt. 1) then
C         call mexErrMsgIdAndTxt ('MATLAB:convhp:nOutput',
C     +                           'One output required.')
C      endif

C     Validate inputs
      mx  = mxGetM(prhs(1))
      nx  = mxGetN(prhs(1))
      my  = mxGetM(prhs(2))
      ny  = mxGetN(prhs(2))
      nkn = mxGetN(prhs(3))

C     Check that inputs are row vectors.
C      if(mx .ne. 1 .or. my .ne. 1) then
C         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector',
C     +                           'Inputs must be row vectors.')
C     Check size of the inputs.
C      elseif(nx .gt. 100 .or. ny .gt. 100) then
C         call mexErrMsgIdAndTxt ('MATLAB:convec:MaxElements',
C     +                 'Inputs must have less than 100 elements.')
C     Check if inputs are numeric.
C      elseif ((mxIsNumeric(prhs(1)) .ne. 1) .or.
C     +        (mxIsNumeric(prhs(2)) .ne. 1)) then
C         call mexErrMsgIdAndTxt ('MATLAB:convec:NonNumeric',
C     +                           'Inputs must be numeric.')
C      endif

C     Create the output array.
      one = 1
      plhs(1)     = mxCreateDoubleMatrix(one,9, 0)
      x_pr        = mxGetPr(prhs(1))
      y_pr        = mxGetPr(prhs(2))
      kn_pr       = mxGetPr(prhs(3))
      k_pr        = mxGetPr(prhs(4))
      j_pr        = mxGetPr(prhs(5))
      nd_pr       = mxGetPr(prhs(6))
      z_pr        = mxGetPr(plhs(1))
      
C     Load the data into Fortran arrays.
      call mxCopyPtrToReal8(x_pr,x,nx*mx)
      call mxCopyPtrToReal8(y_pr,y,ny*my)
      call mxCopyPtrToReal8(kn_pr,kn,nkn)
      call mxCopyPtrToReal8(k_pr,k,one)
      call mxCopyPtrToReal8(j_pr,j,one)
      call mxCopyPtrToReal8(nd_pr,nd,one)

C     Call the computational subroutine.
      call sintpoly3(x,y,kn,k,j,nd,z,nx,ny,mx,my,nkn)

C     Load the output into a MATLAB array.
      call mxCopyReal8ToPtr(z,z_pr,9)

      return
      end
