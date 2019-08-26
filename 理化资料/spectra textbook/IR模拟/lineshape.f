      program main
        implicit none
        integer, parameter:: N=25000
        integer, parameter:: jmax=10000  !step
        real*8:: m_x(0:N-1)
        real*8:: m_y(0:N-1)
        real*8:: m_z(0:N-1)   ! ang
        real*8:: omega(0:N-1)  !cm-1
        real*8:: T1(0:N-1)   !second
        real*8 :: rephi(0:jmax-1)
        real*8 :: imphi(0:jmax-1)
        real*8 :: phi(2*jmax)
        integer::io
        integer::i
        io=25 
        open(unit=io,file="data.csv")
        call readin(N,io,m_x,m_y,m_z,omega,T1)
        omega=omega+1767.2d0
        omega=omega*0.02997924579999719d12
        T1=1.0d-12 
        ! omega=1700.0d0
        ! do i=0,N-1
        ! write(*,*)m_x(i),m_y(i),m_z(i),omega(i)
        ! enddo
        call calcphi(N,10,jmax,m_x,m_y,m_z,omega,T1,rephi,imphi)
        ! do i=0,jmax-1
        ! write(*,*)i, rephi(i), imphi(i)
        ! enddo

        do i=1,jmax
        phi(2*i-1)=rephi(i-1)
        phi(2*i)=imphi(i-1)
        enddo

        call fourier(rephi,imphi,jmax)

        stop
      end program main

      subroutine readin(N,io,m_x,m_y,m_z,omega,T1)
        implicit none
        integer:: N
        real*8:: m_x(0:N-1)
        real*8:: m_y(0:N-1)
        real*8:: m_z(0:N-1)
        real*8:: omega(0:N-1)
        real*8:: T1(0:N-1)
        integer::io
        integer:: i
        rewind(io)
        do i=0,N-1
        read(io,*)m_x(i),m_y(i),m_z(i),omega(i),T1(i)
        call normalize(m_x(i),m_y(i),m_z(i))
        enddo
        return
      end subroutine readin

      subroutine calcphi(N,M,jmax,m_x,m_y,m_z,omega,T1,rephi,imphi)
        implicit none
        integer::N,M,jmax
        ! N: length of the MD trajetory
        ! M: length of the block
        ! jmax: lenght of phi
        real*8 :: dt=1.d-15
        real*8:: m_x(0:N-1)
        real*8:: m_y(0:N-1)
        real*8:: m_z(0:N-1)
        real*8::T1(0:N-1)
        real*8::omega(0:N-1)
        real*8 :: rephi(0:jmax-1)
        real*8 :: imphi(0:jmax-1)
        real*8:: intomega,intT1
        common /dt/ dt

        integer::i,j,k
        integer::ii

        open(unit=77, file="correlation-1ps-rot.dat")
        
        do j=0,jmax-1
        k=0
        rephi(j)=0.0d0
        imphi(j)=0.0d0
        
        do i=0,(N-j-1)/M ! essemble average
        k=k+1

        ! integral here we use some accelarated algorithm
        if (j.le.M .or. i.eq.0) then
          ! calculate the integral normally
          intomega=0.0d0
          intT1=0.0d0
          do ii=i*M,i*M+j
          intomega=intomega+dt*omega(ii)
          intT1=intT1+dt/2.0d0/T1(ii)
          enddo
        else
          ! We keep the integral from the last step
          ! get the overlapped part
          do ii=(i-1)*M,i*M-1
          intomega=intomega-dt*omega(ii)
          intT1=intT1-dt/2.0d0/T1(ii)
          enddo
          ! add the later part
          do ii=(i-1)*M+j+1,i*M+j
          intomega=intomega+dt*omega(ii)
          intT1=intT1+dt/2.0d0/T1(ii)
          enddo
        endif

        ! Euler's formula: exp(ix)=cos(x)+i*sin(x)
         rephi(j)=rephi(j)
     &  +(m_x(i*M+j)*m_x(i*M)+m_x(i*M+j)*m_x(i*M)+m_z(i*M+j)*m_z(i*M))
     &  *dexp(-intT1)*dcos(intomega)
        imphi(j)=imphi(j)
     &  -(m_x(i*M+j)*m_x(i*M)+m_x(i*M+j)*m_x(i*M)+m_z(i*M+j)*m_z(i*M))
     &  *dexp(-intT1)*dsin(intomega)

        enddo
        rephi(j)=rephi(j)/k
        imphi(j)=imphi(j)/k
        if(mod(j,100).eq.0) write(*,*)j,rephi(j),imphi(j)
        write(77,*)j,rephi(j),imphi(j)

        enddo

        return
      end subroutine calcphi

      subroutine fourier(rephi,imphi,nn)
        implicit none
        integer:: nn
        real*8 :: rephi(nn),imphi(nn)
        real*8 :: omega, dt, re,im, omegak
        common /dt/ dt
        integer :: i
        open(unit=99, file="lineshape-1ps-rot.dat")
        ! do omega=-270.0d0,150.0d0,0.01d0
        do omega= 1400.0d0,2000.0d0,1.0d0    !!!!! range
        omegak=omega*0.02997924579999719d12
        re=0.0d0
        im=0.0d0
        do i=1,nn
        re=re
     &  +dcos(omegak*i*dt)*rephi(i)
     &  -dsin(omegak*i*dt)*imphi(i)
        im=im
     &  +dcos(omegak*i*dt)*imphi(i)
     &  +dsin(omegak*i*dt)*rephi(i)
        enddo
        re=re/dsqrt(nn*1.0d0)
        im=im/dsqrt(nn*1.0d0)
        write(99,*) omega, re, im
        enddo
        return
      end subroutine


      subroutine normalize(x,y,z)
        implicit none
        real*8 :: x,y,z,w
        w=x**2+y**2+z**2
        w=dsqrt(w)
        x=x/w
        y=y/w
        z=z/w
        return
      end subroutine normalize
