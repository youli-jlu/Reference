      include 'pesh2o.f'
      include 'parmar.f'
**********************************************************************
      implicit real*8(a-h,o-z)
      parameter (klmax=2**20,mv=9)
      character(24) system1
      character(24) system2
      dimension trad1(nr1max,nr1max),xrad1(nr1max)
      dimension trad2(nr2max,nr2max),xrad2(nr2max)
      dimension wvib1(nr2max,nr2max),evib1(nr2max)
      dimension wvib2(nr2max,nr2max),evib2(nr2max)
      dimension phi0a(nr1max,nr2max)
      dimension wave(nr1max,nr2max),wave0(mv,nr1max,nr2max)
      dimension wz(klmax,mv),w(mv),cw(mv)
      common /bdata/ama,amb,amc,rbohr,const,pi
      common /phi/phi0(nr1max,nr2max),phi1(nr1max,nr2max),
     & vp(nr1max,nr2max)
       ma=max0(nr1max,nr2max)
      data rbohr/0.5291772108d0/
      data ama,amb,amc,rbohr/1.00782503207d0,15.99491461956d0,
     & 1.00782503207d0, 0.5291772108d0/
      open(10,file='vector.dat')     ! the eigenvector in the Lanczos representation
      open(20,file='value.dat')      ! the eigenvalue of vector
      open(40,file='compare.dat')    ! save the compare values between the
                                     ! true eigenvalue and the calculated value
      open(50,file='wave.dat')       ! the wavefunction with r1 and r2 
      open(71,file='r1.dat')         ! save the radial grid points of r1
      open(72,file='r2.dat')         ! save the radial grid points of r2
      open(74,file='ayz.dat')
      open(6,file='podvr-11.out')
      open(8,file='vdw.inp')
      open(30,file='ck1.dat',form='unformatted')
      const=2.0d0*109737.31568525d0
c      call potread
      nr1=70
      nr2=70
      write(*,*)'nr2,nr1=',nr2,nr1
      amau=1822.887427d0
      ama=ama*amau
      amb=amb*amau
      amc=amc*amau
      amu1=1.0d0*amau
      amu2=1.0d0*amau
      amu3=1.0d0*amau
c       damu1=1/ama+1/(amb+amc)
c       damu2=1/amb+1/amc
c       amu1=1/damu1
c       amu2=1/damu2
      pi=dacos(-1.0d0)
c      write(*,*)'JTOT=',jtot,'  PARITY=',iparity
c     kl : the total number of Lanczos steps
      rmin1=-1.0
      rmax1=1.0
      rmin2=-1.0
      rmax2=1.0
      npd=200
      jtot=0
      iparity=0
      if (jtot.lt.iparity) stop 'You make a big mistake'
      kl=250
      nv=3
      vcut=5.0/27.212d0
      r1e=0.0d0
      r2e=0.0d0
      write(*,*)'r1e,r2e= ',r1e,r2e
c     read the eigenvectors from file
      do i=1,kl
      do j=1,nv
         read(10,*) wz(i,j)
      enddo
      enddo
c  read the eigenvalue from file
      do i=1,nv
      read(20,*)k,w(k)
      enddo
c  end read

c  set the initial wavefunction
      do m=1,nv
      do i=1,nr1
      do j=1,nr2
         wave0(m,i,j)=0.0d0
      enddo
      enddo
      enddo

c      call kinet(nr1,r1min,r1max,nr1max,trad1,xrad1,amu1)
c      write(*,*) xrad1(1)
c      write(*,*) xrad1(nr1)

      call podvr1(npd,rmin1,rmax1,r2e,wvib1,evib1,trad1,xrad1,
     & nr1,nr2,amu1)

      call podvr2(npd,rmin2,rmax2,r1e,wvib2,evib2,trad2,xrad2,
     & nr1,nr2,amu2)

      vmin=100000.0
        do 3 i=1,nr1
           r1=xrad1(i)
          do 3 j=1,nr2
            r2=xrad2(j)                         
            vpot=v(r1,r2)
            if (vpot.gt.vcut) vpot=vcut
            vp(i,j)=vpot
            if (vpot.lt.vmin)then
           vmin=vpot
           a1=r1
           a2=r2
           endif
  3   continue
      write(6,*) 'vmin=',vmin
      write(*,*)  a1,a2
      vmin1=vmin*const
      write(*,*) 'vmin(cm-1)=',vmin1
      write(*,*) 'initial state'
      idum=-156356
      sum=0.0
      do 10 i=1,nr1
      do 10 j=1,nr2
            x=ran2(idum)
            ww=2*x-1
            phi0(i,j)=ww
            sum=sum+ww**2
10    continue
      sum=dsqrt(sum)
      do 20 i=1,nr1
      do 20 j=1,nr2
            phi0(i,j)=phi0(i,j)/sum
  20  continue

c     Lanczos propagation
      write(*,*)' propagation begin'
      beta=0.0
      nmtot=nr1*nr2
      do 100 ll=1,kl
      do 103 m=1,nv
      do 104 i=1,nr1
      do 104 j=1,nr2
           wave0(m,i,j)=wave0(m,i,j)+phi0(i,j)*wz(ll,m)
104       continue
103    continue

c     pw1 is the overlap between  the Lanczos state and the matrix which
c    elements are all equal to 1.
      pw1=0
      do 11 i=1,nr1
      do 11 j=1,nr2
       pw1=pw1+phi0(i,j)
 11    continue
c     calculate H*Phi

      call hphi(xrad1,xrad2,trad1,trad2,nr1,nr2,amu1,amu2,ma)

c     calculate the coefficient alpha
       alpha=0
      do 12 i=1,nr1
      do 12 j=1,nr2
      alpha=alpha+phi0(i,j)*(phi1(i,j)-beta*phi0a(i,j))
 12    continue

c     calculate coefficient beta
      beta1=0
      do 13 i=1,nr1
      do 13 j=1,nr2
      phi1(i,j)=phi1(i,j)-alpha*phi0(i,j)-beta*phi0a(i,j)
      beta1=beta1+phi1(i,j)**2
13    continue
      beta=dsqrt(beta1)
      do 14 i=1,nr1
      do 14 j=1,nr2
       phi0a(i,j)=phi0(i,j)
       phi0(i,j)=phi1(i,j)/beta
  14   continue
      write(6,1001)ll,alpha,beta,pw1
      write(30)ll,alpha,beta,pw1
1001  format(1x,i6,3f20.13)
 100  continue
      do 200 m=1,nv
           s=0.0d0
           do 201 i=1,nr1
           do 201 j=1,nr2
               wave(i,j)=wave0(m,i,j)
               s=s+wave(i,j)**2
 201      continue
         write(50+m,*) 'VARIABLES = "Q1/bohr", "Q3/bohr", "w"   '
         write(50+m,*) 'ZONE I=70, J=70, F=point'
           do i=1,nr1
           do j=1,nr2
                wave(i,j)=wave(i,j)/dsqrt(s)
           write(50+m,*)xrad1(i),xrad2(j),wave(i,j) 
           end do
           end do
         call part4(nmtot,wave,phi0)
         call hphi(xrad1,xrad2,trad1,trad2,nr1,nr2,amu1,amu2,ma)
         call part5(nmtot,phi1,wave,w0)
         cw(m)=w0*const
 200     continue
         do 2 m=1,nv
            write(40,1002) m,w(m),cw(m),cw(m)-w(m)
  2     continue
1000    format(1x,40f9.6)
1002    format(1x,i5,3f18.7)
      end 

      subroutine kinet(n,a,b,md,t,x,am)
      implicit real*8(a-h,o-z)
      dimension t(md,md),x(md)
      pi=dacos(-1.0d0)
      n1=n+1
      dh=(b-a)/n1
      x(1)=a+dh
      do 1 i=2,n
1      x(i)=x(i-1)+dh
      do 2 i=1,n
           w=(2*n1**2+1)/3.0d0-1/dsin(pi*i/n1)**2
           t(i,i)=w*pi**2/(4*am*(b-a)**2)
      do 2 j=1,i-1
           w=1/dsin(pi*(i-j)/(2*n1))**2-1/dsin(pi*(i+j)/(2*n1))**2
           t(i,j)=w*pi**2*(-1)**(i-j)/(4*am*(b-a)**2)
           t(j,i)=t(i,j)
2     continue
      end

      SUBROUTINE TQL(MD,N,Z,D,E)
C       Z(-1) A  Z  =D                                    
C       A = Z
C       EIGENVALUE         D(I)
C       EIGENFUNCTION      Z(J,I),J=1,N
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION   D(MD),E(MD),Z(MD,MD)
      EPS=1D-12
      NITER=50
      CALL TRED2(MD,N,Z,D,E)
      DO 10 I=2,N
  10  E(I-1)=E(I)
      F=0.0D0
      B=0.0D0
      E(N)=0.0D0
      DO 20 L=1,N
      J=0
      H=EPS*(DABS(D(L))+DABS(E(L)))
      LP1=L+1
      IF (B-H) 30,40,40
  30  B=H
  40  DO 50 M=L,N
      IF (DABS(E(M))-B) 60,60,50
  50  CONTINUE
  60  IF (M-L) 70,80,70
  70  IF (J-NITER) 90,100,90
  90  J=J+1
      P=(D(LP1)-D(L))/(2*E(L))
      R=DSQRT(P*P+1)
      IF (P) 110,111,111
  110 H=D(L)-E(L)/(P-R)
      GOTO 130
  111 H=D(L)-E(L)/(P+R)
  130 DO 140 I=L,N
  140 D(I)=D(I)-H
      F=F+H
      P=D(M)
      C=1.0D0
      S=0.0D0
      MM1=M-1
      IF (MM1-L) 270,280,280
  280 DO 120 LMIP=L,MM1
      I=L+MM1-LMIP
      IP1=I+1
      G=C*E(I)
      H=C*P
      IF (DABS(P)-DABS(E(I))) 160,170,170
  170 C=E(I)/P
      R=DSQRT(C*C+1.0D0)
      E(IP1)=S*P*R
      S=C/R
      C=1.0D0/R
      GOTO 180
  160 C=P/E(I)
      R=DSQRT(C*C+1)
      E(IP1)=S*E(I)*R
      S=1/R                  
      C=C/R
  180 P=C*D(I)-S*G
      D(IP1)=H+S*(C*G+S*D(I))
      DO 190 K=1,N
      H=Z(K,IP1)
      Z(K,IP1)=S*Z(K,I)+C*H
  190 Z(K,I)=C*Z(K,I)-S*H
  120 CONTINUE
  270 E(L)=S*P
      D(L)=C*P
      IF (DABS(E(L))-B) 80,80,70
  80  D(L)=D(L)+F
  20  CONTINUE
      DO 112 I=1,N
      IP1=I+1
      K=I
      P=D(I)
      IF (N-I) 230,230,300
  300 DO 210 J=IP1,N
      IF (D(J)-P) 220,210,210
  220 K=J
      P=D(J)
  210 CONTINUE
  230 IF (K-I) 240,112,240
  240 D(K)=D(I)
      D(I)=P
      DO 260 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
  260 Z(J,K)=P
  112 CONTINUE
      RETURN
  100 STOP '  FAIL'
      END


      SUBROUTINE TRED2(MD,N,Z,D,E)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  D(MD),E(MD),Z(MD,MD)
      BETA=1D-20
      DO 20 NMIP2=2,N
      I=N+2-NMIP2
      IM1=I-1
      IM2=I-2
      L=IM2
      F=Z(I,IM1)
      G=0.0D0
      IF (L) 30,30,40
  40  DO 50 K=1,L
  50  G=G+Z(I,K)*Z(I,K)
  30  H=G+F*F
      IF (G-BETA) 60,60,70
  60  E(I)=F
      H=0.0D0
      GOTO 180
  70  L=L+1
      IF (F) 80,90,90
  90  E(I)=-DSQRT(H)
      G=E(I)
      GOTO 100
  80  E(I)=DSQRT(H)
      G=E(I)
 100  H=H-F*G
      Z(I,IM1)=F-G
      F=0.0D0
      DO 110 J=1,L
      Z(J,I)=Z(I,J)/H
      G=0.0D0
      DO 201 K=1,J
  201 G=G+Z(J,K)*Z(I,K)
      JP1=J+1
      IF (JP1-L) 130,130,140
  130 DO 120 K=JP1,L
  120 G=G+Z(K,J)*Z(I,K)
  140 E(J)=G/H
      F=F+G*Z(J,I)
  110 CONTINUE
      HH=F/(H+H)
      DO 160    J=1,L
      F=Z(I,J)
      E(J)=E(J)-HH*F
      G=E(J)
      DO 170 K=1,J
  170 Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  160 CONTINUE
  180 D(I)=H
  20  CONTINUE
      D(1)=0.0D0
      E(1)=0.0D0
      DO 190 I=1,N
      L=I-1
      IF (D(I)) 202,210,202
  202 IF (L) 210,210,220
  220 DO 230 J=1,L
      G=0.0D0
      DO 240 K=1,L
  240 G=G+Z(I,K)*Z(K,J)
      DO 250 K=1,L
  250 Z(K,J)=Z(K,J)-G*Z(K,I)
  230 CONTINUE
  210 D(I)=Z(I,I)
      Z(I,I)=1.0D0
      IF (L) 260,260,270
  270 DO 280 J=1,L
      Z(I,J)=0.0D0
  280 Z(J,I)=0.0D0
  260 CONTINUE
  190 CONTINUE
      return
      END      

c****************************************************************************
c
      FUNCTION ran2(idum)
c
c----------------------------------------------------------------------------
c   Long period (>2*10^18) random number generator of L'Ecuyer with Bays-
c   Durham shuffle and added safeguards. Returns a uniform random deviate
c   Between 0.0 and 1.0 (exclusive of the endpoints values 0 and 1). Call
c   with idum a negative integer to initialize; thereafter, do not alter
c   idum between successive devistes in a sequence. RNMX should approximate
c   the largest floating value that is less than one.
c****************************************************************************
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
      idum=max(-idum,1)
      idum2=idum
      do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
11    continue
      iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=dmin1(AM*iy,RNMX)
      return
      END


      subroutine sum1(n,x,s1)
      real*8 x(*),s1
      s1=0
      do 1 i=1,n
1     s1=s1+x(i)
      end


      subroutine sum2(n,x1,x2,x3,a,b)
      real*8 x1(*),x2(*),x3(*),a,b
      a=0
      do 1 i=1,n
1     a=a+x2(i)*(x1(i)-b*x3(i))
      end



      subroutine  sum3(n,x1,x2,x3,a,b,b1)
      real*8 x1(*),x2(*),x3(*),a,b,b1
      b1=0
      do 1 i=1,n
           x1(i)=x1(i)-a*x2(i)-b*x3(i)
           b1=b1+x1(i)**2
1     continue
      end


      subroutine part1(n,x1,x2,x3,c)
      real*8 x1(*),x2(*),x3(*),c
      do 1 i=1,n
           x3(i)=x2(i)
           x2(i)=x1(i)/c
1     continue
      end


      subroutine part2(n,x1,x2,x3)
      real*8 x1(*),x2(*),x3(*)
      do 1 i=1,n
1     x1(i)=x3(i)*x2(i)
      end

      subroutine part3(n,m,i,x1,x2,x3,k)
      real*8 x1(m,*),x2(*),x3(k,m)
      do 1 j=1,n
 1    x1(i,j)=x1(i,j)+x2(j)*x3(k,i)
      end
      
      subroutine part4(n,x1,x2)
      real*8 x1(*),x2(*)
      do 1 j=1,n
      x2(j)=x1(j)
 1    continue
      end
         
      subroutine part5(n,x1,x2,x3)
      real*8 x1(*),x2(*),x3
      x3=0.0d0
      do 1 i=1,n 
      x3=x3+x1(i)*x2(i)
 1    continue
        end


      subroutine  hphi(xrad1,xrad2,trad1,trad2,nr1,nr2,amu1,amu2,ma)
      include 'parmar.f'
      implicit real*8(a-h,o-z)
      parameter (nmax=300)
      dimension ww(nmax)
      dimension xrad1(nr1max),trad1(nr1max,nr1max)
      dimension xrad2(nr2max),trad2(nr2max,nr2max)
      common /bdata/ama,amb,amc,rbohr,const,pi
      common /phi/phi0(nr1max,nr2max),phi1(nr1max,nr2max),
     & vp(nr1max,nr2max)

      do 2 i=1,nr1
      do 2 j=1,nr2
   2  phi1(i,j)=vp(i,j)*phi0(i,j)

      do 20 ir2=1,nr2
      do 21 ir1=1,nr1
  21   ww(ir1)=phi0(ir1,ir2)
      do 20 ir1=1,nr1
      sum=0
      do 22 jr1=1,nr1
  22  sum=sum+trad1(ir1,jr1)*ww(jr1)
  20  phi1(ir1,ir2)=phi1(ir1,ir2)+sum

      do 30 ir1=1,nr1
      do 31 ir2=1,nr2
  31  ww(ir2)=phi0(ir1,ir2)
      do 30 ir2=1,nr2
      sum=0
      do 32 jr2=1,nr2
  32  sum=sum+trad2(ir2,jr2)*ww(jr2)
  30  phi1(ir1,ir2)=phi1(ir1,ir2)+sum
       end



      function v(r1,r2)
      implicit real*8(a-h,o-z)
      pi=dacos(-1.0d0)
      Q1=r1
      Q2=0.0
      Q3=r2
      call POTH2O(vpot,Q1,Q2,Q3)
      v=vpot
c     write(*,*) r1,r2,cthi,v
      end

      subroutine podvr1(npd,r1a,r1b,r2e,wvib1,evib1,
     & trad1,xrad1,nr1,nr2,amu1)
      include 'parmar.f'
      implicit real*8(a-h,o-z)
      parameter (npod=300)
c   Calculate PODVR points and kinetic matrix
c   Parameter npod is a local varaiable only used in this subroutine,
c    which value must greater than npd
c   Parameter nr2 must be the same value as in other subroutines.
c   input: npd = the number of primitive sine-DVR points (e.g. npd=150),
c          nr2 = the number of podvr points (eg. nr1=80)
c        [r2a,r2b] the inteval of r1 of interest,
c         r1e = the equilibrium value of r2 
c         cthe = the equilibrium value of cos(theta) in Radau coordinates
c  output: 
c        wvib is the value of the contracted basis at podvr points,
c        evib is the eigenvalue of the contracted basis  ,
c        trad2 is the kinetic matrix in PODVR
c        xrad2 is the PODVR points
c        amu2 is the reduced mass in the kinetic operator related to r2
       dimension evib1(nr1max),wvib1(nr1max,nr1max),wsr(nr1max)
       dimension trad1(nr1max,nr1max),xrad1(nr1max),www(nr1max)
       dimension ht1(nr1max,nr1max)
       dimension trad(npod,npod),xrad(npod),ww(npod),vv(npod)
       dimension trad0(npod,npod)
       dimension t2(npod,npod),dd(npod),cc(npod,npod),tt(npod,npod)
       common /bdata/ama,amb,amc,rbohr,const,pi
       if (npd.gt.npod) stop ' error: npod must greater than npd'
       if (npd.lt.nr2) stop ' erro: nr1 must less than npd'
c       original dvr point
        call kinet(npd,r1a,r1b,npod,trad,xrad,amu1)
c       calculate eigenvalue and eigenvector
        do 1 ir1=1,npd
        do 1 jr1=1,npd
         trad0(ir1,jr1)=trad(ir1,jr1)
 1      t2(ir1,jr1)=trad(ir1,jr1)
c       TT:  sinc function
        do 3 ir1=1,npd
        sum=0
        do 4 jr1=1,npd
        x=dsin(ir1*pi*(xrad(jr1)-r1a)/(r1b-r1a))
        sum=sum+x*x
 4      tt(ir1,jr1)=x
        sum=dsqrt(sum)
        do 3 jr1=1,npd
 3      tt(ir1,jr1)=tt(ir1,jr1)/sum
        
        do 10 ir1=1,npd
        r1=xrad(ir1)
        v1=v(r1,r2e)
        t2(ir1,ir1)=t2(ir1,ir1)+v1
        vv(ir1)=v1
        write(6,*)ir1,r1,v1*const
 10     continue
        call tql(npod,npd,t2,dd,ww)
          write(6,*)' dd',(dd(i)*const,i=1,10)
        do 5 i=1,nr2
  5     evib1(i)=dd(i)
         do 8 i=1,npd
         do 8 j=1,npd
         sum=0
         do 9 k=1,npd
  9      sum=sum+tt(i,k)*t2(k,j)
  8      cc(i,j)=sum
        
        do 11 i=1,nr1
        do 11 j=1,nr1
        sum0=0
        sum=0 
        do 12 k=1,npd
        sum0=sum0+t2(k,i)*vv(k)*t2(k,j)
 12     sum=sum+t2(k,i)*xrad(k)*t2(k,j)
        trad(i,j)=-sum0
 11     ht1(i,j)=sum
        do 13 i=1,nr1
 13     trad(i,i)=trad(i,i)+dd(i)
        call tql(nr1max,nr1,ht1,xrad1,www)
        write(6,*)' xrad2',(xrad1(i),i=1,nr2)
         do 18 n=1,nr1
         sumt=0
         do 38 k=1,nr1
         r1=xrad1(k)
         sum=0 
         do 19 j=1,npd
         term=dsin(j*pi*(r1-r1a)/(r1b-r1a))
 19      sum=sum+cc(j,n)*term
         sumt=sumt+sum**2
 38      wvib1(n,k)=sum
         sumt=dsqrt(sumt)
         do 18 k=1,nr1
 18      wvib1(n,k)=wvib1(n,k)/sumt
         do 24 k=1,nr1
         wss=ht1(1,k)/wvib1(1,k)
         if (wss.lt.0) then
          wss=-wss
          do 29 n=1,nr1
 29       ht1(n,k)=-ht1(n,k)
          endif
          wsr(k)=wss
 24      continue
        do 15 i=1,nr1
        do 15 j=1,nr1
        sum=0 
        do 16 k=1,nr1
  16    sum=sum+trad(i,k)*ht1(k,j)
  15    t2(i,j)=sum
        do 17 i=1,nr1
        do 17 j=1,i
        sum=0 
        do 28 k=1,nr1
  28    sum=sum+ht1(k,i)*t2(k,j)
        trad1(i,j)=sum
        trad1(j,i)=sum
 17      continue
         end
      

      subroutine podvr2(npd,r2a,r2b,r1e,wvib2,evib2,
     & trad2,xrad2,nr1,nr2,amu2)
      include 'parmar.f'
      implicit real*8(a-h,o-z)
      parameter (npod=300)
c   Calculate PODVR points and kinetic matrix
c   Parameter npod is a local varaiable only used in this subroutine,
c    which value must greater than npd
c   Parameter nr2 must be the same value as in other subroutines.
c   input: npd = the number of primitive sine-DVR points (e.g. npd=150),
c          nr2 = the number of podvr points (eg. nr1=80)
c        [r2a,r2b] the inteval of r1 of interest,
c         r1e = the equilibrium value of r2 
c         cthe = the equilibrium value of cos(theta) in Radau coordinates
c  output:
c        wvib is the value of the contracted basis at podvr points,
c        evib is the eigenvalue of the contracted basis  ,
c        trad2 is the kinetic matrix in PODVR
c        xrad2 is the PODVR points
c        amu2 is the reduced mass in the kinetic operator related to r2
       dimension evib2(nr2max),wvib2(nr2max,nr2max),wsr(nr2max)
       dimension trad2(nr2max,nr2max),xrad2(nr2max),www(nr2max)
       dimension ht1(nr2max,nr2max)
       dimension trad(npod,npod),xrad(npod),ww(npod),vv(npod)
       dimension trad0(npod,npod)
       dimension t2(npod,npod),dd(npod),cc(npod,npod),tt(npod,npod)
       common /bdata/ama,amb,amc,rbohr,const,pi
       if (npd.gt.npod) stop ' error: npod must greater than npd'
       if (npd.lt.nr2) stop ' erro: nr1 must less than npd'
c       original dvr point
        call kinet(npd,r2a,r2b,npod,trad,xrad,amu2)
c       calculate eigenvalue and eigenvector
        do 1 ir2=1,npd
        do 1 jr2=1,npd
         trad0(ir2,jr2)=trad(ir2,jr2)
 1      t2(ir2,jr2)=trad(ir2,jr2)
c       TT:  sinc function
        do 3 ir2=1,npd
        sum=0
        do 4 jr2=1,npd
        x=dsin(ir2*pi*(xrad(jr2)-r2a)/(r2b-r2a))
        sum=sum+x*x
 4      tt(ir2,jr2)=x
        sum=dsqrt(sum)
        do 3 jr2=1,npd
 3      tt(ir2,jr2)=tt(ir2,jr2)/sum

        do 10 ir2=1,npd
        r2=xrad(ir2)
        v1=v(r1e,r2) 
        t2(ir2,ir2)=t2(ir2,ir2)+v1
        vv(ir2)=v1
        write(6,*)ir2,r2,v1*const
 10     continue
        call tql(npod,npd,t2,dd,ww)
          write(6,*)' dd',(dd(i)*const,i=1,10)
        do 5 i=1,nr2
  5     evib2(i)=dd(i)
         do 8 i=1,npd
         do 8 j=1,npd
         sum=0
         do 9 k=1,npd
  9      sum=sum+tt(i,k)*t2(k,j)
  8      cc(i,j)=sum

        do 11 i=1,nr2
        do 11 j=1,nr2
        sum0=0
        sum=0
        do 12 k=1,npd
        sum0=sum0+t2(k,i)*vv(k)*t2(k,j)
 12     sum=sum+t2(k,i)*xrad(k)*t2(k,j)
        trad(i,j)=-sum0
 11     ht1(i,j)=sum
        do 13 i=1,nr2
 13     trad(i,i)=trad(i,i)+dd(i)
        call tql(nr2max,nr2,ht1,xrad2,www)
        write(6,*)' xrad2',(xrad2(i),i=1,nr2)
         do 18 n=1,nr2
         sumt=0
         do 38 k=1,nr2
         r2=xrad2(k)
         sum=0
         do 19 j=1,npd
         term=dsin(j*pi*(r2-r2a)/(r2b-r2a))
 19      sum=sum+cc(j,n)*term
         sumt=sumt+sum**2
 38      wvib2(n,k)=sum
         sumt=dsqrt(sumt)
         do 18 k=1,nr2
 18      wvib2(n,k)=wvib2(n,k)/sumt
         do 24 k=1,nr2
         wss=ht1(1,k)/wvib2(1,k)
         if (wss.lt.0) then
          wss=-wss
          do 29 n=1,nr2
 29       ht1(n,k)=-ht1(n,k)
          endif
          wsr(k)=wss
 24      continue
        do 15 i=1,nr2
        do 15 j=1,nr2
        sum=0
        do 16 k=1,nr2
  16    sum=sum+trad(i,k)*ht1(k,j)
  15    t2(i,j)=sum
        do 17 i=1,nr2
        do 17 j=1,i
        sum=0
        do 28 k=1,nr2
  28    sum=sum+ht1(k,i)*t2(k,j)
        trad2(i,j)=sum
        trad2(j,i)=sum
 17      continue
         end






      


