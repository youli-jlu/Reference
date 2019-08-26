C    Calculate the vibrational energy levels of ABA molecule
C    Using Lanczos algorithm and PODVR

C    part 2 of PODVR

C  input: file ' ck.dat'  should be created
C  output: file 'podvr.out'  which containing data :
C       i, de(i),z1(i),z2(3),z3(i)
c    where i is the number of eigenvalues
c      de(i) is the eigenval of the ith root minus the first root
c      z1(i) is the first component of lanczos eigenvector
c      z2(i) is the last component of lanczos eigenvector. This value should
c       be very small (<1e-7) for the converged eigenstates.
c      z3(i) is the overlap between  the eigenstate and a special matrix
c        in which all the elements are equal to 1. For the symmetrical state,
c        this value is large;
c        for the antisymmetrical state, this value is very small (<1e-7 )


	implicit real*8 (a-h,o-z)
	parameter (klmax=2**17, novlp=3)
	dimension diag(klmax), subdiag(klmax), z(klmax,novlp), ztmp(novlp)
       open(25,file='podvr-2.out')
       open(30,file='ck.dat',form='unformatted')
       const=2.0d0*109737.32d0
	write(*,*)'  kl (the Lanczos steps) =?'
!       kl=5000
c	write(*,*)'  kl (the Lanczos steps)=  ',kl
        read(*,*)kl
	write(25,*)'  kl=',kl
	do 1 i=1,kl
 1 	read(30,end=100)k,diag(k),subdiag(k+1),z(k,3)
	matz=1
	matx=novlp-2
	call rst1nx(klmax,kl,diag,subdiag,matz,matx,z,ierr)
	write(25,*) 'kl,ierr=',kl,ierr
	do 2 i=1,kl
 2	DIAG(I)=DIAG(I)*CONST
	do  5 i=1,kl
	write(25,'(i7,f20.12,10e15.6)') i,diag(i),(z(i,j),j=1,2+matx)
  5	continue
	stop
 100	continue
	stop ' not enough data'

	end

      subroutine rst1nx(nm,n,w,e,matz,matx,z,ierr)
c
c	if matz=1, calculate the first, last component of each eigenvector.
c	if matx>0, calculate the overlap between matx vectors and eigenvectors
      integer i,j,n,nm,ierr,matz
      double precision w(n),e(n),z(nm,*)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric tridiagonal matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix.
c
c        w  contains the diagonal elements of the real
c        symmetric tridiagonal matrix.
c
c        e  contains the subdiagonal elements of the matrix in
c        its last n-1 positions.  e(1) is arbitrary.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for imtql1
c           and imtql2.  the normal completion code is zero.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  imtql1(n,w,e,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
crq         do 30 j = 1, n
crq            z(j,i) = 0.0d0
crq   30    continue
c
crq         z(i,i) = 1.0d0
         z(i,1)=0.0d0
         z(i,2)=0.0d0
   40 continue
      z(1,1)=1.d0; z(n,2)=1.d0
c
      call  imtql2(matx,nm,n,w,e,z,ierr)
   50 return
      end
      subroutine imtql1(n,d,e,ierr)
c
      integer i,j,l,m,n,ii,mml,ierr
      double precision d(n),e(n)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure imtql1,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the implicit ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      e(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
c     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = dabs(d(m)) + dabs(d(m+1))
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 215
crq         if (j .eq. 30) go to 1000
         if (j .eq. 30) go to 290
         j = j + 1
c     .......... form shift ..........
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = pythag(g,1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r,g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            r = pythag(f,g)
            e(i+1) = r
            if (r .eq. 0.0d0) go to 210
            s = f / r
            c = g / r
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
c     .......... recover from underflow ..........
  210    d(i+1) = d(i+1) - p
         e(m) = 0.0d0
         go to 105
c     .......... order eigenvalues ..........
  215    if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine imtql2(matx,nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,nm,mml,ierr
crq      double precision d(n),e(n),z(nm,n)
      double precision d(n),e(n),z(nm,*)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure imtql2,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the implicit ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
c     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = dabs(d(m)) + dabs(d(m+1))
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
crq         if (j .eq. 30) go to 1000
         if (j .eq. 30) go to 240
         j = j + 1
c     .......... form shift ..........
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = pythag(g,1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r,g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            r = pythag(f,g)
            e(i+1) = r
            if (r .eq. 0.0d0) go to 210
            s = f / r
            c = g / r
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     .......... form vector ..........
crq         do 180 k = 1, n
crq            f = z(k,i+1)
crq            z(k,i+1) = s * z(k,i) + c * f
crq            z(k,i) = c * z(k,i) - s * f
crq 180     continue
            do 180 k = 1, 2+matx
               f = z(i+1,k)
               z(i+1,k) = s * z(i,k) + c * f
               z(i,k) = c * z(i,k) - s * f
  180       continue
c
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
c     .......... recover from underflow ..........
  210    d(i+1) = d(i+1) - p
         e(m) = 0.0d0
         go to 105
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
c        do 280 j = 1, n
c           p = z(j,i)
c           z(j,i) = z(j,k)
c           z(j,k) = p
c 280    continue
         do 280 j = 1, 2+matx
            p = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end




