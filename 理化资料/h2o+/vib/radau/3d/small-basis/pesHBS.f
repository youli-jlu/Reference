c   ab inito potential energy surface for the ground electronic state of n2o 
c      Please add the following statement in your main program:
c       call potread
c     use the following statement to calculate the PES when needed:
c       call HBSpes(rr1,rr2,theta,v)
c      where rr1 is R(N-N) and rr2 is the R(N-O) in angstron,
c       theta is the enclosed angle in degree
c       v is the potential in hartree

c     implicit real*8(a-h,o-z)
c     iv=2
c     call potread(iv)
c     pi=dacos(-1.0d0)
c     rr1=3.0d0
c     rr2=3.0d0
c     theta=180.0d0
c     call hbspes(iv,rr1,rr2,theta,v)
c     write(*,*) rr1,rr2,theta,v
c     end


     
      subroutine potread(iv)
      implicit real*8(a-h,o-z)
      parameter (pi=3.141592653589793d0)    
      parameter (nr1=24,nr2=25,nth=20)
      common /bpescut/vmin,vasy,vcut
      common /pes/pesmin,ri(nr1),rj(nr2),thek(nth),
     &vcid(3,nr1,nr2,nth)
      data vmin/-423.56309393d0/
      data vasy/-422.50593275d0/
      open(98,file='hbs-pes-acvtz-f12.dat')
c      open(99,file='test.dat')
      vcut=6.0/27.21138d0
      read(98,*)
      pesmin=10000
      do k=1,10
       kk=0
      do i=1,nr1
      do j=1,nr2
        read(98,*)thek(k),ri(i),rj(j),(vcid(m,i,j,k),m=1,3)
        if (vcid(iv,i,j,k).lt.pesmin) then
         pesmin=vcid(iv,i,j,k)
         r1m=ri(i)
         r2m=rj(j)
         r3m=thek(k)
         endif
         if ((vcid(iv,i,j,k)-vmin).gt.vcut) vcid(iv,i,j,k)=vmin+vcut
         kk=21-k
         thek(kk)=360.0d0-thek(k)
         vcid(iv,i,j,kk)=vcid(iv,i,j,k)
        enddo;enddo;enddo
        write(*,*)r1m,r2m,r3m,pesmin
11      continue

        do k=1,nth
        do i=1,nr1
        do j=1,nr2
           rr1=ri(i)
           rr2=rj(j)
           thth=thek(k)
         vcid(iv,i,j,k)=vcid(iv,i,j,k)-vmin
c        write(99,1001)thth,rr1,rr2,vcid(iv,i,j,k)
c 1001   format (1x,2f15.6,f8.4,f15.8)
      enddo;enddo;enddo
      end



      subroutine hbspes(iv,r10,r20,theta,v)
      implicit real*8(a-h,o-z)
      parameter (nr1=24,nr2=25,nth=20)    
      common /bpescut/vmin,vasy,vcut
      common /pes/pesmin,ri(nr1),rj(nr2),thek(nth),
     &vcid(3,nr1,nr2,nth)
      parameter (rbohr=0.5291772108d0)
      parameter (pi=3.141592653589793d0)
      rr1=r10
      rr2=r20
      thth=theta
      if (rr1.lt.1.2d0) rr1=1.2d0
      if (rr2.lt.2.1d0) rr2=2.1d0
      if (rr1.gt.7.0d0) rr1=7.0d0
      if (rr2.gt.8.0d0) rr2=8.0d0
      if (thth.lt.90.0) thth=90.0
      call SPl3(iv,rr1,rr2,thth,vg)
      v=vg
      end

      subroutine spl3(iv,r10,r20,th,vg)
      implicit real*8(a-h,o-z)
      parameter (nr1=24,nr2=25,nth=20)    
      parameter (m=100,n=100,l=100,ir2=2)
       dimension xt(m)  
       dimension dty(l),ddty(l),s1(1),ds1(1),dds1(1),h1(l)
       dimension dny(n),ddny(n),s2(1),ds2(1),dds2(1),h2(n)
       dimension dhy(m),ddhy(m),s3(1),ds3(1),dds3(1),h3(m)
       dimension y(m),ss(m),sss(m),y2(l)
      common /pes/pesmin,ri(nr1),rj(nr2),thek(nth),
     & vcid(3,nr1,nr2,nth)
      data dy1,dyn/1.0d30,1.0d30/
      do 20 i=1,nr1
      do 10 j=1,nr2
         nh=0
      do 2 k=1,nth
         nh=nh+1
         xt(nh)=thek(k)
         y(nh)=vcid(iv,i,j,k)
   2  continue
        if (th.lt.xt(1)) th=xt(1)
      call spline(xt,y,nh,dy1,dyn,y2)
      call splint(xt,y,y2,nh,th,y3)
 22   continue
        ss(j)=y3
   10   continue
   
      nr=0
      do 5 j=1,nr2
        nr=nr+1
        xt(nr)=rj(j)
        y(nr)=ss(j)
        y(nr)=(xt(nr)**ir2)*y(nr) 
   5  continue
      call spline(xt,y,nr,dy1,dyn,y2)
      call splint(xt,y,y2,nr,r20,yw2)
 33    continue
       sss(i)=yw2/r20**ir2
   20    continue
  
        nr=0
        do i=1,nr1
        nr=nr+1
        xt(nr)=ri(i)
        y(nr)=sss(i)
        y(nr)=(xt(nr)**ir2)*y(nr) 
        enddo
        call spline(xt,y,nr,dy1,dyn,y2)
        call splint(xt,y,y2,nr,r10,yw3)
   44   continue
        vg=yw3/r10**ir2
        end


          SUBROUTINE ESPL2(X,Y,N,DY1,DYN,XX,M,DY,DDY,S,DS,DDS,T,H)
           parameter (nmax=100)
!           DIMENSION X(N),Y(N),XX(M),DY(N),DDY(N)
!        DIMENSION S(M),DS(M),DDS(M),H(N)
           dimension x(nmax),y(nmax),xx(nmax),dy(nmax),ddy(nmax)
           dimension s(nmax),ds(nmax),dds(nmax),h(nmax)
	   DOUBLE PRECISION X,Y,XX,DY,DDY,S,DS,DDS,H,DY1,DYN,
     *                   T,H0,H1,BETA,ALPHA
c	   DY(1)=-0.5
	   H0=X(2)-X(1)
c	   H(1)=3.0*(Y(2)-Y(1))/(2.0*H0)-DY1*H0/4.0
		dy(1)=0
		h(1)=0 
	   DO 10 J=2,N-1
	  H1=X(J+1)-X(J)
	  ALPHA=H0/(H0+H1)
	  BETA=(1.0-ALPHA)*(Y(J)-Y(J-1))/H0
	  BETA=3.0*(BETA+ALPHA*(Y(J+1)-Y(J))/H1)
	  DY(J)=-ALPHA/(2.0+(1.0-ALPHA)*DY(J-1))
	  H(J)=(BETA-(1.0-ALPHA)*H(J-1))
	  H(J)=H(J)/(2.0+(1.0-ALPHA)*DY(J-1))
	  H0=H1
10	  CONTINUE
c	  DY(N)=(3.0*(Y(N)-Y(N-1))/H1+DYN*H1/2.0-H(N-1))
c     *         /(2.0+DY(N-1))
	  DY(N)=0
	  DO 20 J=N-1,1,-1
20	  DY(J)=DY(J)*DY(J+1)+H(J)
	  DO 30 J=1,N-1
30	  H(J)=X(J+1)-X(J)
	  DO 40 J=1,N-1
	  H1=H(J)*H(J)
	  DDY(J)=6.0*(Y(J+1)-Y(J))/H1-
     *           2.0*(2.0*DY(J)+DY(J+1))/H(J)
40	  CONTINUE
	  H1=H(N-1)*H(N-1)
	  DDY(N)=6.0*(Y(N-1)-Y(N))/H1+
     *            2.0*(2.0*DY(N)+DY(N-1))/H(N-1)
	  T=0.0
	  DO 50 I=1,N-1
	  H1=0.5*H(I)*(Y(I)+Y(I+1))
	  H1=H1-H(I)*H(I)*H(I)*(DDY(I)+DDY(I+1))/24.0
	  T=T+H1
50	  CONTINUE
	  DO 70 J=1,M
	  IF (XX(J).GE.X(N)) THEN
	    I=N-1
	  ELSE
	    I=1
60	    IF (XX(J).GT.X(I+1)) THEN
	      I=I+1
	      GOTO 60
	    END IF
	  END IF
	  H1=(X(I+1)-XX(J))/H(I)
	  S(J)=(3.0*H1*H1-2.0*H1*H1*H1)*Y(I)
	  S(J)=S(J)+H(I)*(H1*H1-H1*H1*H1)*DY(I)
	  DS(J)=6.0*(H1*H1-H1)*Y(I)/H(I)
	  DS(J)=DS(J)+(3.0*H1*H1-2.0*H1)*DY(I)
	  DDS(J)=(6.0-12.0*H1)*Y(I)/(H(I)*H(I))
	  DDS(J)=DDS(J)+(2.0-6.0*H1)*DY(I)/H(I)
	  H1=(XX(J)-X(I))/H(I)
	  S(J)=S(J)+(3.0*H1*H1-2.0*H1*H1*H1)*Y(I+1)
	  S(J)=S(J)-H(I)*(H1*H1-H1*H1*H1)*DY(I+1)
	  DS(J)=DS(J)-6.0*(H1*H1-H1)*Y(I+1)/H(I)
	  DS(J)=DS(J)+(3.0*H1*H1-2.0*H1)*DY(I+1)
	  DDS(J)=DDS(J)+(6.0-12.0*H1)*Y(I+1)/(H(I)*H(I))
	  DDS(J)=DDS(J)-(2.0-6.0*H1)*DY(I+1)/H(I)
70	  CONTINUE
	  RETURN
	  END


C##################################################################
C# SPLINE ROUTINES
C#            Numerical recipes in fortran
C#            Cambrige University Press
C#            York, 2nd edition, 1992.
C##################################################################
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      implicit double precision  (a-h,o-z)
      DIMENSION xa(n),y2a(n),ya(n)
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.0d0) write(6,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.0d0
      return
      END
C##############################################################################
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit double precision  (a-h,o-z)
      DIMENSION x(n),y(n),y2(n)
      PARAMETER (NMAX=100)
      DIMENSION u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0d0
        y2(i)=(sig-1.0d0)/p
        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.0d0
        un=0.0d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

