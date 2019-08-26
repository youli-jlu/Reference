      include 'pesh2o.f'
      include 'parmar.f'
**********************************************************************
!This program is used to calculate the ro-vibrational energy levels of
!HBS system
!everytime I seperate the calculation to two parts,odd and even parity
!for a fixed J.
c    Calculate the vibrational energy levels of ABC molecule
C    Using Lanczos algorithm and PODVR
C    part 1 of PODVR
c    all input data must be in Radau coordinates.
c    so all the data in bl-ba coordinates will be changed
c    to Radau coordinates.
c    the radial input data are in Radau coordinats.but the angle
c    is in bl-ba coordinates.
c    nr1 is the points of H-B,nr2 is the points of B-S.
      implicit real*8(a-h,o-z)
c     character(24) system1
c     character(24) system2
      parameter (klmax=2**17,nth=120,nr1=70,nr2=70,mv=5)
c    nth,nr1,nr2 must be the actural used values
      dimension trad1(nr1max,nr1max),xrad1(nr1max)
      dimension wvib1(nr1max,nr1max),evib1(nr1max)
      dimension trad2(nr2max,nr2max),xrad2(nr2max)
      dimension wvib2(nr2max,nr2max),evib2(nr2max)
      dimension tang(0:jmax,nthmax,nthmax),glx(0:jmax,nthmax)
      dimension ht1(0:jmax,nthmax,nthmax)
      dimension ht2(0:jmax,nthmax,nthmax)
      dimension phi0a(0:jmax,nr1max,nr2max,nthmax) 
      dimension wave(0:jmax,nr1max,nr2max,nthmax)
      dimension wave0(mv,0:jmax,nr1,nr2,nth)
      dimension wz(klmax,mv),w(mv),cw(mv)
      common /bdata/ama,amb,amc,rbohr,const,pi
      common /phi/phi0(0:jmax,nr1max,nr2max,nthmax),
     &phi1(0:jmax,nr1max,nr2max,nthmax),
     & vp(0:jmax,nr1max,nr2max,nthmax)
      data ama,amb,amc,rbohr/1.00782503207d0,15.99491461956d0,
     & 1.00782503207d0, 0.5291772108d0/
c      system1=ctime(time())
c      print*,' The program began at  ',system1
      open(10,file='vector.dat')     ! the eigenvector in the Lanczos representation
      open(20,file='value.dat')      ! the eigenvalue of vector
      open(40,file='compare.dat')    ! save the compare values between the
                                     ! true eigenvalue and the calculated value
      open(50,file='wavek.dat')      ! the wavefunction fixed the angle
      open(60,file='wavei.dat')      ! the wavefunction fixed the r2
      open(70,file='wavej.dat')      ! the wavefunction fixed the r1
      open(71,file='r1.dat')         ! save the radial grid points of r1 
      open(72,file='r2.dat')         ! save the radial grid points of r2 
      open(73,file='angle.dat')      ! save the angle grid points
      open(74,file='ayz.dat')
      open(6,file='podvr-11.out')
      open(30,file='ck1.dat',form='unformatted')
c      call potread(2)
      const=2.0d0*109737.32d0
      amau=1822.887427d0
      ama=ama*amau
      amb=amb*amau
      amc=amc*amau
      amu1=ama
      amu2=amc
c	damu1=1/ama+1/(amb+amc)
c	damu2=1/amb+1/amc
c	amu1=1/damu1
c	amu2=1/damu2
      pi=dacos(-1.0d0)
c     the calculation parameters
c     npd: the sine DVR points for R1
c     kl : the total number of Lanczos steps
        
      rmin1=1.0
      rmax1=7.5
      rmin2=1.0
      rmax2=7.5
      npd=200
      jtot=0
      iparity=0
      if (jtot.lt.iparity) stop 'You make a big mistake'
      kl=400
      nv=5
      vcut=10.0/27.212d0
      rr1=1.8102073945d0
      rr2=1.8102073945d0
      theta=104.5d0
      cthae=dcos(theta/180.0d0)
      write(*,*) 'R(H-O)= ',rr1
      write(*,*) 'R(O-H)= ',rr2
      write(*,*) 'A(H-O-H)= ',theta
      call bondtorad(rr1,rr2,cthae,r1e,r2e,cthe)
      write(*,*)'r1e,r2e,cthe= ',r1e,r2e,cthe
c     call kinet(nr1,rmin1,rmax1,nr1max,trad1,xrad1,amu1)
      call podvr1(npd,rmin1,rmax1,r2e,cthe,wvib1,evib1,
     &                  trad1,xrad1,amu1)
      do i=1,nr1
      write(71,*) xrad1(i)
      enddo
      call podvr2(npd,rmin2,rmax2,r1e,cthe,wvib2,evib2,
     &trad2,xrad2,amu2)
      do i=1,nr2
      write(72,*) xrad2(i)
      enddo
c     call kinet(nr2,rmin2,rmax2,nr2,trad2,xrad2,amc)
      call angheg(glx,tang,ht1,ht2,jtot,iparity,nth)
      write(*,*)'angle part ends'
      do k=1,nth
        if (mod(jmod,2).eq.0) then
            write(73,*) glx(0,k)
            else
            write(73,*) glx(1,k)
            endif
      enddo
c     read the eigenvectors from file
      do i=1,kl
         do j=1,nv
             read(10,*) wz(i,j)
             enddo
             enddo
c     read the eigenvalue from file
      do i=1,nv
          read(20,*)k,w(k)
          enddo
c     calculate potential matrix   
      vmin=100000
      jmod=jtot+iparity
      if (mod(jmod,2).eq.0) then
      do 3 mz=0,jtot
      do 3 i=1,nr1
      r1=xrad1(i)
      tjot=(jtot*(jtot+1)-2*mz**2)/(2*amu1*r1**2)
      do 3 j=1,nr2
      r2=xrad2(j)
      do 3 k=1,nth
       cthi=glx(mz,k)
      call radtobond(r1,r2,cthi,rr1,rr2,cth)
        vpot=v(rr1,rr2,cth)
        if (vpot.gt.vcut) vpot=vcut
        vp(mz,i,j,k)=vpot+tjot
        cthi=dacos(cth)*180.d0/pi
        if (vpot.lt.vmin) then
        vmin=vpot
        a1=r1
        a2=r2
        a3=cthi
        kbend=k
        irad=i
        jrad=j
      endif
   3  continue
      write(6,*)'  vmin=',vmin
      write(*,*)'  vmin=',vmin
      write(6,*)  a1,a2,a3,i,j,k
c     initial state  
      idum=-156356
      sum=0.0
      do 10 i=1,nr1
      do 10 j=1,nr2
      do 10 k=1,nth
      do 10 mz=0,jtot
      x=ran2(idum)
      ww=2*x-1
      phi0(mz,i,j,k)=ww
      sum=sum+ww**2
  10  continue
      sum=dsqrt(sum)
      do 20 i=1,nr1
      do 20 j=1,nr2
      do 20 k=1,nth
      do 20 mz=0,jtot
        phi0(mz,i,j,k)=phi0(mz,i,j,k)/sum
 20   continue

c     Lanczos  propagation
      write(*,*)' propagation begin'
      beta=0.0d0
      nmtot=nr1*nr2*nth*(jtot+1)
      do 100 ll=1,kl

       do 110 l=1,nv
         do 111 k=1,nth
         do 111 j=1,nr2
         do 111 i=1,nr1
         do 111 mz=0,jtot
           wave0(l,mz,i,j,k)=wave0(l,mz,i,j,k)+phi0(mz,i,j,k)*wz(ll,l)
 111       continue
 110      continue
c     pw1 is the overlap between  the Lanczos state and the matrix which 
c     elements are all equal to 1.
      pw1=0.0d0
      do 11 i=1,nr1
      do 11 j=1,nr2
      do 11 k=1,nth
      do 11 mz=0,jtot
       pw1=pw1+phi0(mz,i,j,k)
 11   continue
c     intead of calculating 'pw1' above,
c     we use the subroutine ' sum1'  to save CPU time.
c     call sum1(nmtot,phi0,pw1)
c     calculate H*Phi
       call hphi(xrad1,xrad2,trad1,trad2,tang,ht1,ht2,iparity,jtot,
     & amu1,amu2)
c     call hphi(xrad1,xrad2,trad1,trad2,tang)
c     calculate the coefficient alpha
      alpha=0.0d0
      do 12 i=1,nr1
      do 12 j=1,nr2
      do 12 k=1,nth
      do 12 mz=0,jtot
       alpha=alpha+phi0(mz,i,j,k)*(phi1(mz,i,j,k)-beta*phi0a(mz,i,j,k))
 12   continue
c	call sum2(nmtot,phi1,phi0,phi0a,alpha,beta)
c     calculate coefficient beta
      beta1=0.0d0
      do 13 i=1,nr1
      do 13 j=1,nr2
      do 13 k=1,nth
      do 13 mz=0,jtot
       phi1(mz,i,j,k)=phi1(mz,i,j,k)-alpha*phi0(mz,i,j,k)-
     &beta*phi0a(mz,i,j,k)
      beta1=beta1+phi1(mz,i,j,k)**2
  13  continue
c	call sum3(nmtot,phi1,phi0,phi0a,alpha,beta,beta1)
      beta=dsqrt(beta1)
      do 14 i=1,nr1
      do 14 j=1,nr2
      do 14 k=1,nth 
      do 14 mz=0,jtot
      phi0a(mz,i,j,k)=phi0(mz,i,j,k)
      phi0(mz,i,j,k)=phi1(mz,i,j,k)/beta
  14  continue
c	call part1(nmtot,phi1,phi0,phi0a,beta) 
      write(6,1001)ll,alpha,beta,pw1
      write(30)ll,alpha,beta,pw1
  100 continue
      do 300 l=1,nv
         s=0.0d0
         do 301 i=1,nr1
         do 301 j=1,nr2
         do 301 k=1,nth
         do 301 mz=0,jtot
           wave(mz,i,j,k)=wave0(l,mz,i,j,k)
          s=s+wave(mz,i,j,k)**2
 301     continue
         do 302 i=1,nr1
         do 302 j=1,nr2
         do 302 k=1,nth
         do 302 mz=0,jtot
            wave(mz,i,j,k)=wave(mz,i,j,k)/dsqrt(s)
 302     continue
         kki=50+l
         write(kki,*) 'VARIABLES = "R1 /a_0", "R2 /a_0", "v/a.u."'
         write(kki,*) 'ZONE I=70, J=70, F=point'
         do 303 i=1,nr1  
         do 303 j=1,nr2  
         write(kki,1004) xrad1(i),xrad2(j),wave(0,i,j,kbend)
 303     continue 

         kkj=60+l
         write(kkj,*) 'VARIABLES = "R1 /a_0", "theta/degree", "v/a.u."'
         write(kkj,*) 'ZONE I=120, J=70, F=point'
         do 304 i=1,nr1
         do 304 k=1,nth
         write(kkj,1004) xrad1(i),dacos(glx(0,k))*180.d0/pi,
     &    wave(0,i,jrad,k)/((1.-glx(0,k)**2)**0.25)
 304     continue

         kkk=80+l
         write(kkk,*) 'VARIABLES = "R2 /a_0", "theta/degree", "v/a.u."'
         write(kkk,*) 'ZONE I=120, J=70, F=point'
         do 305 j=1,nr2
         do 305 k=1,nth 
         write(kkk,1004) xrad2(j),dacos(glx(0,k))*180.d0/pi,
     &    wave(0,irad,j,k)/((1.-glx(0,k)**2)**0.25)
 305     continue
      call part4(nmtot,wave,phi0)
      call hphi(xrad1,xrad2,trad1,trad2,tang,ht1,ht2,iparity,jtot,
     & amu1,amu2)
      call part5(nmtot,phi1,wave,w0)
      cw(l)=w0*const
 300  continue
      do 2 i=1,nv
      write(40,1002) i,w(i),cw(i),cw(i)-w(i)
 2    continue

      else
      do 4  mz=1,jtot
      do 4 i=1,nr1
      r1=xrad1(i)
      tjot=(jtot*(jtot+1)-2*mz**2)/(2*amu1*r1**2)
      do 4 j=1,nr2
      r2=xrad2(j)
      do 4 k=1,nth
       cthi=glx(mz,k)
      call radtobond(r1,r2,cthi,rr1,rr2,cth)
      vpot=v(rr1,rr2,cth)
        if (vpot.gt.vcut) vpot=vcut
        vp(mz,i,j,k)=vpot+tjot
        cthi=dacos(cth)*180.d0/pi
        if (vpot.lt.vmin) then
        vmin=vpot
        a1=r1
        a2=r2
        a3=cthi
      endif
   4   continue
      write(6,*)'  vmin=',vmin
      write(*,*)'  vmin=',vmin
      write(*,*)  a1,a2,a3
c     initial state
      idum=-156356
      sum=0.0d0
      do 30 i=1,nr1
      do 30 j=1,nr2
      do 30 k=1,nth
      do 30 mz=1,jtot
      x=ran2(idum)
      ww=2*x-1
      phi0(mz,i,j,k)=ww
        sum=sum+ww**2
  30  continue
      sum=dsqrt(sum)
      do 40 i=1,nr1
      do 40 j=1,nr2
      do 40 k=1,nth
      do 40 mz=1,jtot
      phi0(mz,i,j,k)=phi0(mz,i,j,k)/sum
 40   continue
c     Lanczos  propagation
      write(*,*)' propagation begin'
      beta=0.0d0
c     nmtot=nr1*nr2*nth
      do 200 ll=1,kl
c     pw1 is the overlap between  the Lanczos state and the matrix which
c     elements are all equal to 1.
      pw1=0
      do 18 i=1,nr1
      do 18 j=1,nr2
      do 18 k=1,nth
      do 18 mz=1,jtot
       pw1=pw1+phi0(mz,i,j,k)
 18    continue
c     intead of calculating 'pw1' above,
c     we use the subroutine ' sum1'  to save CPU time.
c     call sum1(nmtot,phi0,pw1)
c     calculate H*Phi
        call hphi(xrad1,xrad2,trad1,trad2,tang,ht1,ht2,iparity,jtot,
     & amu1,amu2)
c     call hphi(xrad1,xrad2,trad1,trad2,tang)
c     calculate the coefficient alpha
      alpha=0.0d0
      do 15 i=1,nr1
      do 15 j=1,nr2
      do 15 k=1,nth
      do 15 mz=1,jtot
       alpha=alpha+phi0(mz,i,j,k)*(phi1(mz,i,j,k)-beta*phi0a(mz,i,j,k))
 15   continue
c     call sum2(nmtot,phi1,phi0,phi0a,alpha,beta)
c     calculate coefficient beta
      beta1=0
      do 16 i=1,nr1
      do 16 j=1,nr2
      do 16 k=1,nth
      do 16 mz=1,jtot
       phi1(mz,i,j,k)=phi1(mz,i,j,k)-alpha*phi0(mz,i,j,k)-
     &beta*phi0a(mz,i,j,k)
        beta1=beta1+phi1(mz,i,j,k)**2
  16  continue
c     call sum3(nmtot,phi1,phi0,phi0a,alpha,beta,beta1)
      beta=dsqrt(beta1)
      do 17 i=1,nr1
      do 17 j=1,nr2
      do 17 k=1,nth
      do 17 mz=1,jtot
      phi0a(mz,i,j,k)=phi0(mz,i,j,k)
      phi0(mz,i,j,k)=phi1(mz,i,j,k)/beta
  17  continue
c     call part1(nmtot,phi1,phi0,phi0a,beta)
      write(6,1001)ll,alpha,beta,pw1
      write(30)ll,alpha,beta,pw1
  200   continue
      do 400 l=1,nv
         s=0.0d0
         do 401 i=1,nr1
         do 401 j=1,nr2
         do 401 k=1,nth
         do 401 mz=1,jtot
           wave(mz,i,j,k)=wave0(l,mz,i,j,k)
          s=s+wave(mz,i,j,k)**2
 401     continue
         do 402 i=1,nr1
         do 402 j=1,nr2
         do 402 k=1,nth
         do 402 mz=1,jtot
            wave(mz,i,j,k)=wave(mz,i,j,k)/dsqrt(s)
 402     continue
         kki=50+l
         write(kki,*) 'VARIABLES = "R1 /a_0", "R2 /a_0", "v/a.u."'
         write(kki,*) 'ZONE I=70, J=70, F=point'
         do 403 i=1,nr1
         do 403 j=1,nr2
         write(kki,1004) xrad1(i),xrad2(j),
     &   wave(1,i,j,kbend)
 403     continue

         kki=60+l
         write(kki,*) 'VARIABLES = "R1 /a_0", "theta/degree", "v/a.u."'
         write(kki,*) 'ZONE I=120, J=70, F=point'
         do 404 i=1,nr1
         do 404 k=1,nth
         write(kki,1004) xrad1(i),dacos(glx(1,k))*180.d0/pi,
     &    wave(1,i,jrad,k)/((1.-glx(1,k)**2)**0.25)
 404     continue
         kki=80+l
         write(kki,*) 'VARIABLES = "R2 /a_0", "theta/degree", "v/a.u."'
         write(kki,*) 'ZONE I=120, J=70, F=point'
         do 405 j=1,nr2
         do 405 k=1,nth
         write(kki,1004) xrad2(j),dacos(glx(1,k))*180.d0/pi,
     &    wave(1,irad,j,k)/((1.-glx(1,k)**2)**0.25) 
 405     continue
      call part4(nmtot,wave,phi0)
      call hphi(xrad1,xrad2,trad1,trad2,tang,ht1,ht2,iparity,jtot,
     & amu1,amu2)
      call part5(nmtot,phi1,wave,w0)
      cw(l)=w0*const
 400  continue
      do 6 i=1,nv
      write(40,1002) i,w(i),cw(i),cw(i)-w(i)
 6    continue
      end if 

1001  format(1x,i5,4f20.13)
1002  format(1x,i5,3f18.7)
1004  format(1x,10f15.7)

c     system2=ctime(time())
c     print*,' The program ended at   ',system2
      end

      subroutine sum1(n,x,s1)
      real*8 x(*),s1
      s1=0
      do 1 i=1,n
 1    s1=s1+x(i)
      end

	subroutine sum2(n,x1,x2,x3,a,b)
	real*8 x1(*),x2(*),x3(*),a,b
	a=0
	do 1 i=1,n
 1	a=a+x2(i)*(x1(i)-b*x3(i))
	end


	subroutine  sum3(n,x1,x2,x3,a,b,b1)
	real*8 x1(*),x2(*),x3(*),a,b,b1
	b1=0
	do 1 i=1,n
	x1(i)=x1(i)-a*x2(i)-b*x3(i)
	b1=b1+x1(i)**2
 1	continue
	end
	
	subroutine part1(n,x1,x2,x3,c)
	real*8 x1(*),x2(*),x3(*),c
	do 1 i=1,n
	x3(i)=x2(i)
	x2(i)=x1(i)/c
 1	continue
	end

	subroutine part2(n,x1,x2,x3)
	real*8 x1(*),x2(*),x3(*)
	do 1 i=1,n
 1	x1(i)=x3(i)*x2(i)
	end


c      subroutine part3(n,m,i,x1,x2,x3,k)
c      real*8 x1(m,*),x2(*),x3(k,m)
c      do 1 j=1,n
c   1   x1(i,j)=x1(i,j)+x2(j)*x3(k,i)
c      end


      subroutine part4(n,x1,x2)
      real*8 x1(*),x2(*)
      do 1 j=1,n
      x2(j)=x1(j)
   1  continue
      end

      subroutine part5(n,x1,x2,x3)
      real*8 x1(*),x2(*),x3
      x3=0.0d0
      do 1 i=1,n
      x3=x3+x1(i)*x2(i)
   1  continue
      end





      subroutine hphi(xrad1,xrad2,trad1,trad2,tang,ht1,ht2,iparity,jtot,
     &amu1,amu2)
c      subroutine hphi(xrad1,xrad2,trad1,trad2,tang)
      include 'parmar.f'
      implicit real*8(a-h,o-z)
c      calculate H*Phi
C      input:
C      nr1,nr2,nth are the numbers of DVR points in r1,r2 and theta
C      xrad(i): the DVR points in r1 and r2 (for ABC molecule)
c      trad(i,j): The DVR kinetic matrix for r1 and r2
c      tang(i,j): The DVR kinetic matrix for theta
C      important notes for use:
C      (a) the parameter values of nth and nr1 must be the same as
C          in the main program;
C      (b) The two common statements must also appeared in the main program,
C          ama,amb,amc are atomic mass (in atomic unit) of atoms A,B,C
C     rbohr=0.5291772108d0, const=2.0d0*109737.32d0,pi=dacos(-1.0d0)
C     Phi1=H*Phi0, vp is the potential matrix

      parameter (nth=120,nr1=70,nr2=70)
      parameter (nmax=300)
      dimension tang(0:jmax,nthmax,nthmax),ww(nthmax)
      dimension xrad1(nr1max),trad1(nr1max,nr1max)
      dimension xrad2(nr2max),trad2(nr2max,nr2max)
      dimension ht1(0:jmax,nthmax,nthmax)
      dimension ht2(0:jmax,nthmax,nthmax)
      common /bdata/ama,amb,amc,rbohr,const,pi
      common /phi/phi0(0:jmax,nr1max,nr2max,nthmax),
     &phi1(0:jmax,nr1max,nr2max,nthmax),
     &    vp(0:jmax,nr1max,nr2max,nthmax)
        
      jmod=jtot+iparity
c	nmtot=nr1*nr2*nth
c	call part2(nmtot,phi1,phi0,vp)
      if(mod(jmod,2).eq.0) then 
      do 10 mz=0,jtot
      do 10 ir1=1,nr1
      do 10 ir2=1,nr2
      do 10 ith=1,nth
  10  phi1(mz,ir1,ir2,ith)=vp(mz,ir1,ir2,ith)*phi0(mz,ir1,ir2,ith)
      do 20 mz=0,jtot
      do 20 ir2=1,nr2
      do 20 ith=1,nth
      do 21 ir1=1,nr1
  21  ww(ir1)=phi0(mz,ir1,ir2,ith)

      do 20 ir1=1,nr1
      sum=0
      do 22 jr1=1,nr1
 22   sum=sum+trad1(ir1,jr1)*ww(jr1)
 20   phi1(mz,ir1,ir2,ith)=phi1(mz,ir1,ir2,ith)+sum
      do 30 mz=0,jtot
      do 30 ir1=1,nr1
      do 30 ith=1,nth
      do 31 ir2=1,nr2
  31  ww(ir2)=phi0(mz,ir1,ir2,ith)
      do 30 ir2=1,nr2
      sum=0
      do 32 jr2=1,nr2
  32  sum=sum+trad2(ir2,jr2)*ww(jr2)
  30  phi1(mz,ir1,ir2,ith)=phi1(mz,ir1,ir2,ith)+sum
      do 40 mz=0,jtot
      do 40 ir1=1,nr1
      r1=xrad1(ir1)
      do 40 ir2=1,nr2
      r2=xrad2(ir2)
      term=1/(2*amu1*r1**2)+1/(2*amu2*r2**2)
      do 41 ith=1,nth
 41   ww(ith)=phi0(mz,ir1,ir2,ith)
      do 40 ith=1,nth
      sum=0
      do 42 jth=1,nth
 42   sum=sum+tang(mz,ith,jth)*ww(jth)
 40   phi1(mz,ir1,ir2,ith)=phi1(mz,ir1,ir2,ith)+sum*term
!calculate h(k+1,k)*phi(k)
      do 70 mz=0,jtot-1
      do 70 ir1=1,nr1
      r1=xrad1(ir1)
      term=-1.0d0/(2.0d0*amu1*r1**2)
      do 70 ir2=1,nr2
      do 71 ith1=1,nth
         ww(ith1)=phi0(mz,ir1,ir2,ith1)
  71  continue
      do 70 ith=1,nth
      sum=0.0d0
      do 72 jth=1,nth
      sum=sum+ht1(mz,ith,jth)*ww(jth)
  72  continue
      phi1(mz+1,ir1,ir2,ith)=phi1(mz+1,ir1,ir2,ith)+sum*term
  70  continue
!calculate h(k-1,k)*phi(k)
      do 80 mz=1,jtot
      do 80 ir1=1,nr1
      r1=xrad1(ir1)
      term=-1.0/(2.0d0*amu1*r1**2)
      do 80 ir2=1,nr2
      do 81 ith1=1,nth
 81   ww(ith1)=phi0(mz,ir1,ir2,ith1)
      do 80 ith=1,nth
      sum=0.0d0
      do 82 jth=1,nth
 82   sum=sum+ht2(mz,ith,jth)*ww(jth)
      phi1(mz-1,ir1,ir2,ith)=phi1(mz-1,ir1,ir2,ith)+sum*term
 80   continue

      else
      do 3  mz=1,jtot
      do 3  ir1=1,nr1
      do 3 ir2=1,nr2
      do 3  ith=1,nth
   3     phi1(mz,ir1,ir2,ith)=vp(mz,ir1,ir2,ith)*phi0(mz,ir1,ir2,ith)
      do 50 mz=1,jtot
      do 50 ir2=1,nr2
      do 50 ith=1,nth
      do 51 ir1=1,nr1
  51  ww(ir1)=phi0(mz,ir1,ir2,ith)
      do 50 ir1=1,nr1
      sum=0
      do 52 jr1=1,nr1
 52   sum=sum+trad1(ir1,jr1)*ww(jr1)
 50   phi1(mz,ir1,ir2,ith)=phi1(mz,ir1,ir2,ith)+sum
      do 60 mz=1,jtot
      do 60 ir1=1,nr1
      do 60 ith=1,nth
      do 61 ir2=1,nr2
  61  ww(ir2)=phi0(mz,ir1,ir2,ith)
      do 60 ir2=1,nr2
      sum=0
      do 62 jr2=1,nr2
  62  sum=sum+trad2(ir2,jr2)*ww(jr2)
  60  phi1(mz,ir1,ir2,ith)=phi1(mz,ir1,ir2,ith)+sum
      do 7 mz=1,jtot
      do 7 ir1=1,nr1
        r1=xrad1(ir1)
      do 7 ir2=1,nr2
        r2=xrad2(ir2)
      term=1/(2*amu1*r1**2)+1/(2*amu2*r2**2)
      do 8 ith=1,nth
   8   ww(ith)=phi0(mz,ir1,ir2,ith)
      do 7 ith=1,nth
      sum=0
      do 9 jth=1,nth
   9   sum=sum+tang(mz,ith,jth)*ww(jth)
   7   phi1(mz,ir1,ir2,ith)=phi1(mz,ir1,ir2,ith)+sum*term
      if (jtot.gt.1) then
!calculate h(k+1,k)*phi(k)
      do 90 mz=1,jtot-1
      do 90 ir1=1,nr1
      r1=xrad1(ir1)
      term=-1.0d0/(2.0d0*amu1*r1**2)
      do 90 ir2=1,nr2
      do 91 ith1=1,nth
         ww(ith1)=phi0(mz,ir1,ir2,ith1)
  91  continue
      do 90 ith=1,nth
      sum=0.0d0
      do 92 jth=1,nth
      sum=sum+ht1(mz,ith,jth)*ww(jth)
  92  continue
      phi1(mz+1,ir1,ir2,ith)=phi1(mz+1,ir1,ir2,ith)+sum*term
  90  continue
!calculate h(k-1,k)*phi(k)
      do 100 mz=2,jtot
      do 100 ir1=1,nr1
      r1=xrad1(ir1)
      term=-1.0/(2.0d0*amu1*r1**2)
      do 100 ir2=1,nr2
      do 101 ith1=1,nth
 101   ww(ith1)=phi0(mz,ir1,ir2,ith1)
      do 100 ith=1,nth
      sum=0.0d0
      do 102 jth=1,nth
 102  sum=sum+ht2(mz,ith,jth)*ww(jth)
      phi1(mz-1,ir1,ir2,ith)=phi1(mz-1,ir1,ir2,ith)+sum*term
 100  continue
      end if     
      end if 
      end

      subroutine kinet(n,a,b,md,t,x,am)
      implicit real*8(a-h,o-z)
      dimension t(md,md),x(md)
      pi=dacos(-1.0d0)
      n1=n+1
      dh=(b-a)/n1
      x(1)=a+dh
      do 1 i=2,n
  1   x(i)=x(i-1)+dh
      do 2 i=1,n
      w=(2*n1**2+1)/3.0d0-1/dsin(pi*i/n1)**2
      t(i,i)=w*pi**2/(4*am*(b-a)**2)
      do 2 j=1,i-1
      w=1/dsin(pi*(i-j)/(2*n1))**2-1/dsin(pi*(i+j)/(2*n1))**2
      t(i,j)=w*pi**2*(-1)**(i-j)/(4*am*(b-a)**2)
      t(j,i)=t(i,j)
 2    continue
      end


      subroutine bondtorad(r1,r2,cth,rr1,rr2,ct)
      implicit real*8(a-h,o-z)
c     transfer from bondlength-bond angle coordiantes to Radau coordinates
c     input:
c         r1 and r2 are in bohr, cth=cos(phi)
c     output:
c     Radau coordinates: rr1 and rr2 in bohr, ct=cos(theta)
      common /bdata/ama,amb,amc,rbohr,const,pi
      alp=dsqrt(amb/(ama+amb+amc))
      am12=ama+amc
      b1=(alp-1)*ama/am12+1
      b2=(alp-1)*amc/am12
      rr1=dsqrt(b1**2*r1**2+b2**2*r2**2+2*b1*b2*r1*r2*cth)
      rr2=dsqrt((b1-1)**2*r1**2+
     &    (b2+1)**2*r2**2+2*(b1-1)*(b2+1)*r1*r2*cth)
      ct=(b1*(b1-1)*r1**2+b2*(b2+1)*r2**2+(b1*(1+b2)+
     &  b2*(b1-1))*r1*r2*cth)/(rr1*rr2)
      end

      subroutine radtobond(r1,r2,cth,rr1,rr2,ct)
      implicit real*8(a-h,o-z)
c     transfer from Radaui coordinate to bondlength-bond angle coordiantes i
c     input:
c     Radau coordinate:  r1 and r2 are in bohr, cth=cos(theta)
c     output:
c     rr1 and rr2 are in bohr, ct=cos(phi)
      common /bdata/ama,amb,amc,rbohr,const,pi
      alp=dsqrt(amb/(ama+amb+amc))
      am12=ama+amc
      b1=(1/alp-1)*ama/am12+1
      b2=(1/alp-1)*amc/am12
      rr1=dsqrt(b1**2*r1**2+b2**2*r2**2+2*b1*b2*r1*r2*cth)
      rr2=dsqrt((b1-1)**2*r1**2+
     &    (b2+1)**2*r2**2+2*(b1-1)*(b2+1)*r1*r2*cth)
      ct=(b1*(b1-1)*r1**2+b2*(b2+1)*r2**2+(b1*(1+b2)+
     &  b2*(b1-1))*r1*r2*cth)/(rr1*rr2)
      end


      subroutine angheg(glx,tang,ht1,ht2,jtot,iparity,nth)
      include 'parmar.f'
      implicit real*8(a-h,o-z)
      dimension tang(0:jmax,nthmax,nthmax),glx(0:jmax,nthmax)
      dimension ht1(0:jmax,nthmax,nthmax)
      dimension ht2(0:jmax,nthmax,nthmax)
      dimension ht3(nthmax,nthmax),htt3(0:jmax,0:nthmax,nthmax)
      dimension glx1(nthmax),tang1(nthmax,nthmax)
      jmod=jtot+iparity
        if (mod(jmod,2).eq.0) then
         do mz=0,jtot
          nth1=nth       
            call gsleg(nthmax,nth1,glx1,tang1,mz,ht3)
            do j=1,nth
               glx(mz,j)=glx1(j)
                   do i=1,nth
                      ii=mz+i-1
                      tang(mz,i,j)=tang1(i,j)
                      htt3(mz,ii,j)=ht3(i,j)
                   end do
            end do
          end do

      do 4 mz=0,jtot-1
      if (mz.eq.0) then
      term=dsqrt(2.0d0*(jtot*(jtot+1.0d0)-mz*(mz+1.0d0)))
      else
      term=dsqrt(jtot*(jtot+1.0d0)-mz*(mz+1.0d0))
      end if
      do 4 ith1=1,nth
      do 4 ith2=1,nth
      sum=0
      do 5 j=mz+1,nth+mz-1    
      a=dsqrt(j*(j+1.0d0)-mz*(mz+1.0d0))
      sum=sum+a*htt3(mz+1,j,ith1)*htt3(mz,j,ith2)
  5   continue
      ht1(mz,ith1,ith2)=sum*term
  4   continue
      do 6 mz=1,jtot
      if (mz.eq.1) then
      term=dsqrt(2.0d0*(jtot*(jtot+1.0d0)-mz*(mz-1.0d0)))
      else
      term=dsqrt(jtot*(jtot+1.0d0)-mz*(mz-1.0d0))
      end if
      do 6 ith1=1,nth
      do 6 ith2=1,nth
        sum=0.0
      do 7 j=mz,nth+mz-2
              a=dsqrt(j*(j+1.0d0)-mz*(mz-1.0d0))
  7   sum=sum+a*htt3(mz-1,j,ith1)*htt3(mz,j,ith2)
      ht2(mz,ith1,ith2)=sum*term
  6   continue
      else
          do mz=1,jtot
          nth1=nth
           call gsleg(nthmax,nth1,glx1,tang1,mz,ht3)
            do j=1,nth
               glx(mz,j)=glx1(j)
                   do i=1,nth
                      ii=mz+i-1
                      tang(mz,i,j)=tang1(i,j)
                      htt3(mz,ii,j)=ht3(i,j)
                   end do
            end do
          end do
      do 8 mz=1,jtot-1
      if (mz.eq.0) then
      term=dsqrt(2.0d0*(jtot*(jtot+1.0d0)-mz*(mz+1.0d0)))
      else
      term=dsqrt(jtot*(jtot+1.0d0)-mz*(mz+1.0d0))
      end if
      do 8 ith1=1,nth
       do 8 ith2=1,nth
      sum=0
      do 9 j=mz+1,nth+mz-1
      a=dsqrt(j*(j+1.0d0)-mz*(mz+1.0d0))
      sum=sum+a*htt3(mz+1,j,ith1)*htt3(mz,j,ith2)
  9   continue
      ht1(mz,ith1,ith2)=sum*term
  8   continue
      do 10 mz=2,jtot
      if (mz.eq.1) then
      term=dsqrt(2.0d0*(jtot*(jtot+1.0d0)-mz*(mz-1.0d0)))
      else
      term=dsqrt(jtot*(jtot+1.0d0)-mz*(mz-1.0d0))
      end if
      do 10 ith1=1,nth
      do 10 ith2=1,nth
      sum=0.0
      do 11 j=mz,nth+mz-2
      a=dsqrt(j*(j+1.0d0)-mz*(mz-1.0d0))
  11   sum=sum+a*htt3(mz-1,j,ith1)*htt3(mz,j,ith2)
      ht2(mz,ith1,ith2)=sum*term
  10   continue
      end if
      end

 
      SUBROUTINE gsleg(md,nth,glx,t,mz,ht3)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (legm=300)
      dimension D(legm),E(legm),ht3(md,md),pl(0:legm),wr(legm,legm)
      dimension glx(md),t(md,md),ind(legm),ws(md),htk(legm,legm)
      common /bdata/ama,amb,amc,rbohr,const,pi
      DF(L,M)=DSQRT((L+M+1.0D0)*(L-M+1)/((2*L+1)*(2*L+3)))
      x1=-1.0d0
      x2=1.0d0
      NG3=nth      
      if (ng3.gt.legm) stop ' legm too small'
      DO 2 I=1,NG3
      DO 2 J=1,NG3
  2   HTk(I,J)=0
      DO 3 I=1,NG3-1
  3   HTk(I+1,I)=DF(MZ+I-1,MZ)
      CALL TQL(legm,NG3,HTk,D,E)
      do 11 j=1,ng3
        x=d(j)
        call sylm(legm,ng3,mz,x,pl)
          do 11 i=1,ng3
            wr(i,j)=pl(i-1)
  11  continue
      nth=0
      do 5 i=1,ng3
           if (d(i).gt.x1.and.d(i).lt.x2) then
            nth=nth+1
             glx(nth)=d(i)
              ind(nth)=i
             x=d(i)
           if (htk(1,i)/wr(1,i).lt.0) then
            do ii=1,ng3
            htk(ii,i)=-htk(ii,i)
            enddo
           endif
        endif
  5   continue
      if (nth.gt.md) stop ' md too small'
      do 8 k=1,nth
           k1=ind(k)
           ws(k)=htk(1,k1)/wr(1,k1)
      do 8 kk=1,k
           kk1=ind(kk)
           sum=0
            do 9 l=1,ng3
            lz=l+mz-1
9          sum=sum+htk(l,k1)*lz*(lz+1)*htk(l,kk1)
           t(k,kk)=sum
8     t(kk,k)=sum
      do 20 n=1,nth
      do 20 k=1,nth
      ht3(n,k)=htk(n,k)
   20 continue
      end

         subroutine sylm(nd,n,kz,x,ylm)
         implicit real*8(a-h,o-z)
         dimension ylm(0:nd)
         df1(j,m)=dsqrt(((j+1.0d0)**2-kz**2)/((2*j+1.0d0)*(2*j+3.0d0)))  
         if (kz.eq.0) then
         ylm(0)=1/dsqrt(2.0d0)
         else 
         t0=1
         st=dsqrt(1-x**2)
         do 1 m=1,kz
 1       t0=t0*dsqrt((2*m-1.0d0)/(2.0d0*m))
         ylm(0)=(-1)**kz*dsqrt((2*kz+1.0d0)/2.0d0)*t0*st**kz
         endif
         ylm(1)=x*ylm(0)/df1(kz,kz)
         do 2 j=1,n-1
         ylm(j+1)=(x*ylm(j)-df1(j+kz-1,kz)*ylm(j-1))/df1(j+kz,kz)
 2       continue
         end 



      subroutine nlp(md,n,x,pl,mz)
      include 'parmar.f'
      implicit real*8(a-h,o-z)
      dimension pl(0:md),plm(0:md,0:jmax)
      plm(0,0)=1/dsqrt(2.0d0)
      plm(1,0)=x*dsqrt(1.5d0)
      plm(1,1)=dsqrt((1.0d0-x**2)*0.75d0)
      plm(2,1)=x*dsqrt((1-x**2)*15.0d0/4.0d0)
      do 1 i=1,nthmax-1
           a=(i+1)/dsqrt((2*i+1.0d0)*(2*i+3.0d0))
           b=i/dsqrt((2*i-1.0d0)*(2*i+1.0d0))
1     plm(i+1,0)=(x*plm(i,0)-b*plm(i-1,0))/a
	 if(mz.ge.1) then
      do 2 i=2,nthmax-1
      a=dsqrt((i+1+1.0d0)*(i-1+1.0d0)/((2*i+1.0d0)*(2*i+3.0d0)))
      b=dsqrt((i+1.0d0)*(i-1.0d0)/((2*i+1.d0)*(2*i-1.0d0)))
2     plm(i+1,1)=(x*plm(i,1)-b*plm(i-1,1))/a
      if (mz.gt.1) then
          c=dsqrt(1-x**2)
          do i=1,mz-1 
            a=2.0d0*i/dsqrt((mz-i)*(mz+i+1.0d0))
            b=dsqrt((mz+i)*(mz-i+1.0d0)/((mz+i+1.0d0)*(mz-i)))
            plm(mz,i+1)=(a*plm(mz,i)+b*(1-x**2)*plm(mz,i-1))/c
          end do
       do i=1,mz 
       a=2.0d0*i/dsqrt((mz+1.0d0-i)*(mz+i+2.0d0))
       b=dsqrt((mz+i+1.0d0)*(mz-i+2.0d0)/((mz+i+2.0d0)*(mz-i+1.0d0)))
        plm(mz+1,i+1)=(a*plm(mz+1,i)+b*(1-x**2)*plm(mz+1,i-1))/c
       end do
       do i=mz+1,nthmax-1
       a=dsqrt((i+mz+1.0d0)*(i-mz+1.0d0)/((2*i+1.0d0)*(2*i+3.0d0)))
       b=dsqrt((i+mz)*(i-mz)/((2*i+1.0d0)*(2*i-1.0d0)))
       plm(i+1,mz)=(x*plm(i,mz)-b*plm(i-1,mz))/a
       end do
      end if
	end if
      do i=1,n
      pl(i-1)=plm(mz+i-1,mz)
      end do
      end




c****************************************************************************

c

      FUNCTION ran2(idum)

c

c----------------------------------------------------------------------------

c

c   Long period (>2*10^18) random number generator of L'Ecuyer with Bays-

c   Durham shuffle and added safeguards. Returns a uniform random deviate

c   Between 0.0 and 1.0 (exclusive of the endpoints values 0 and 1). Call

c   with idum a negative integer to initialize; thereafter, do not alter

c   idum between successive devistes in a sequence. RNMX should approximate

c   the largest floating value that is less than one.

c

c

c****************************************************************************

c

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

11      continue

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

      function v(r1,r2,cth)
      implicit real*8(a-h,o-z)
      pi=dacos(-1.0d0)
      rr1=r1
      rr2=r2
      theta=cth
      if(theta.gt.1.d0)theta=1.0d0
      if(theta.lt.-1.d0)theta=-1.0d0
      thth=dacos(theta)
      call POTS(vpot,rr1,rr2,thth)
      v=vpot
c     write(*,*) r1,r2,cthi,v
      end

       subroutine podvr1(npd,r1a,r1b,r2e,cthe,wvib,evib,
     &                  trad1,xrad1,amu1)
       include 'parmar.f'
      implicit real*8(a-h,o-z)
      parameter (npod=300)
      parameter (nth=120,nr1=70,nr2=70)
c     Calculate PODVR points and kinetic matrix
c     Parameter npod is a local varaiable only used in this subroutine,
c     which value must greater than npd
c     Parameter nr1 must be the same value as in other subroutines.
c     input:
c     npd = the number of primitive sine-DVR points (e.g. npd=150),
c     nr1 = the number of podvr points (eg. nr1=80)
c     [r1a,r1b] the inteval of r1 of interest,
c     r2e = the equilibrium value of r2 in Radau coordinates
c     cthe = the equilibrium value of cos(theta) in Radau coordinates
c  output: 
c     wvib is the value of the contracted basis at podvr points,
c     evib is the eigenvalue of the contracted basis  ,
c     trad1 is the kinetic matrix in PODVR
c     xrad1 is the PODVR points
c     amu1 is the reduced mass in the kinetic operator related to r1
       dimension evib(nr1max),wvib(nr1max,nr1max),wsr(nr1max)
       dimension trad1(nr1max,nr1max),xrad1(nr1max),www(nr1max)
       dimension ht1(nr1max,nr1max)
       dimension trad(npod,npod),xrad(npod),ww(npod),vv(npod)
       dimension trad0(npod,npod)
       dimension t2(npod,npod),dd(npod),cc(npod,npod),tt(npod,npod)
       common /bdata/ama,amb,amc,rbohr,const,pi
       if (npd.gt.npod) stop ' error: npod must greater than npd'
       if (npd.lt.nr1) stop ' erro: nr1 must less than npd'
c       original dvr point
        call kinet(npd,r1a,r1b,npod,trad,xrad,amu1)
c       calculate eigenvalue and eigenvector
        do 1 ir2=1,npd
        do 1 jr2=1,npd
         trad0(ir2,jr2)=trad(ir2,jr2)
 1      t2(ir2,jr2)=trad(ir2,jr2)
c       TT:  sinc function
        do 3 ir2=1,npd
        sum=0
        do 4 jr2=1,npd
        x=dsin(ir2*pi*(xrad(jr2)-r1a)/(r1b-r1a))
        sum=sum+x*x
 4      tt(ir2,jr2)=x
        sum=dsqrt(sum)
        do 3 jr2=1,npd
 3      tt(ir2,jr2)=tt(ir2,jr2)/sum
 
        do 10 ir2=1,npd
        r1=xrad(ir2)
         call radtobond(r1,r2e,cthe,rr1,rr2,cth)
        v1=v(rr1,rr2,cth)
        t2(ir2,ir2)=t2(ir2,ir2)+v1
        vv(ir2)=v1
        write(6,*)ir2,r1,v1*27.212
 10     continue
        call tql(npod,npd,t2,dd,ww)
          write(6,*)' dd',(dd(i)*const,i=1,10)
        do 5 i=1,nr1
  5     evib(i)=dd(i)
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
        call tql(nr1,nr1,ht1,xrad1,www)
        write(6,*)' xrad1',(xrad1(i),i=1,nr1)
         do 18 n=1,nr1
         sumt=0
         do 38 k=1,nr1
         r2=xrad1(k) 
         sum=0
         do 19 j=1,npd
         term=dsin(j*pi*(r2-r1a)/(r1b-r1a))
 19      sum=sum+cc(j,n)*term
         sumt=sumt+sum**2
 38      wvib(n,k)=sum
         sumt=dsqrt(sumt)
         do 18 k=1,nr1
 18      wvib(n,k)=wvib(n,k)/sumt
         do 24 k=1,nr1
         wss=ht1(1,k)/wvib(1,k)
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



      subroutine podvr2(npd,r2a,r2b,r1e,cthe,wvib2,evib2,
     &trad2,xrad2,amu2)
      include 'parmar.f'
      implicit real*8(a-h,o-z)
      parameter (npod=300)
      parameter (nth=120,nr1=70,nr2=70)
c     Calculate PODVR points and kinetic matrix
c     Parameter npod is a local varaiable only used in this subroutine,
c     which value must greater than npd
c     Parameter nr2 must be the same value as in other subroutines.
c     input:
c     npd = the number of primitive sine-DVR points (e.g. npd=150),
c     nr2 = the number of podvr points (eg. nr1=80)
c     [r2a,r2b] the inteval of r1 of interest,
c     r1e = the equilibrium value of r2 in Radau coordinates
c     cthe = the equilibrium value of cos(theta) in Radau coordinates
c  output: 
c     wvib is the value of the contracted basis at podvr points,
c     evib is the eigenvalue of the contracted basis  ,
c     trad2 is the kinetic matrix in PODVR
c     xrad2 is the PODVR points
c     amu2 is the reduced mass in the kinetic operator related to r2
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
        call radtobond(r1e,r2,cthe,rr1,rr2,cth)
        v1=v(rr1,rr2,cth)
        t2(ir2,ir2)=t2(ir2,ir2)+v1
        vv(ir2)=v1
        write(6,*)ir2,r2,v1*27.212
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
