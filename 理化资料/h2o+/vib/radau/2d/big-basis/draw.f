      
	parameter (nth=120,nr1=70,nr2=nr1)
c      nth,nr1,nr2 must be the actural used values
c	 nr1 is the points of h--xe,nr2 is the points of xe--i.
 
      dimension xrad1(nr1), xrad2(nr2),glx(nth)
      dimension wavek(nr1,nr2),wavei(nr2,nth),wavej(nr1,nth)
c	open(20,file='outdrawth-1.dat')
c	open(30,file='outdrawr2-1.dat')
c	open(40,file='outdrawr1-1.dat')
      open(50,file='wavek.dat')      ! the wavefunction fixed the angle
      open(60,file='wavei.dat')      ! the wavefunction fixed the r2
      open(70,file='wavej.dat')      ! the wavefunction fixed the r1
      open(81,file='r1.dat')         ! save the radial grid points of r1
      open(82,file='r2.dat')         ! save the radial grid points of r2
      open(83,file='angle.dat')     ! save the angle grid points
      pi=dacos(-1.0d0)
      nv=5
	read(81,*) (xrad1(i),i=1,nr1)
	read(82,*) (xrad2(j),j=1,nr2)
	read(83,*) (glx(k),k=1,nth)
      do kk=1,nv
      read(50,1000) ((wavek(i,j),j=1,nr2),i=1,nr1)
      read(60,1004) ((wavei(i,k),i=1,nr1),k=1,nth)
      read(70,1000) ((wavej(j,k),j=1,nr2),k=1,nth)
	

      kkk=20+kk
      write(kkk,*) 'VARIABLES = "r1/a_0", "r2/a_0", "v/a0"'
      write(kkk,*) 'ZONE I=70, J=70, F=point'
      do i=1,nr1
 	do j=1,nr2
 	write(kkk,999)xrad1(i),xrad2(j),wavek(i,j)
    	enddo
	enddo

      kki=30+kk
      write(kki,*) 'VARIABLES = "r1/a_0", "theta", "v/a0"'
      write(kki,*) 'ZONE I=70, J=120, F=point'
      do i=1,nr1
     	do k=1,nth
        th1=acos(glx(k))*180.d0/pi
	write(kki,999)xrad1(i),th1,wavei(i,k)
	enddo
	enddo
      
      kkj=40+kk
      write(kkj,*) 'VARIABLES = "r2/a_0", "theta", "v/a0"'
      write(kkj,*) 'ZONE I=70, J=120, F=point'
      do j=1,nr2
 	do k=1,nth
        th2=acos(glx(k))*180.d0/pi
	write(kkj,999)xrad2(j),th2,wavej(j,k)
	enddo
	enddo
      enddo 
1000    format(1x,40f9.6)
1002    format(1x,i5,3f18.7)
1004    format(1x,60f9.6)
999    format(1x,3f20.8)
      end


