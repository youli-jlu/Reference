!*********************************************************************
!***  This program will calculate the eigenvectors for the selected **
!***  eigenvalues. And it uses the method of inverse iteration.     **
!*********************************************************************
      implicit real*8(a-h,o-z)
      parameter(klmax=2**17,mr=30)
      dimension diag(klmax),subdiag(klmax),w(klmax),wz(klmax,mr)
!open:
!     33) stores the eigenvectors of the selected eigenvalues.
!     34) stores the selected eigenvalues
!     35) stores the unformatted data got from the PODVR-1.f
      open(33,file='vector.dat')
      open(34,file='value.dat')
      open(35,file='ck.dat',form='unformatted')
       const=2.0d0*109737.32d0
      write(*,*)' The lanczos step kl is,kl='
      read(*,*)kl
	write(*,*)' The number of your selected eigenvalues is,n='
      read(*,*)n
c      read(34,*) kl,n
      do 200	i=1,kl
        read(35,end=100)k,diag(k),subdiag(k+1)
200   continue  
      do 201 i=1,n
         read(34,*)k,w(k)
         w(k)=w(k)/const
          
201     continue 
        call tinvit(kl,diag,subdiag,n,w,wz,ierr)
        write(*,*) 'w=',w(1)*const
        if (ierr.ne.0) goto 800
	 
        do i=1,kl
	  do j=1,n
        write(33,*) wz(i,j)
        end do
	enddo
        write(*,*) 'ierr=?',ierr
      stop
 100  continue
      stop ' not enough data'
800   continue
       write(*,*) 'ierr=?',ierr
      stop 'tinvit failed'
      end
    
									

      SUBROUTINE TINVIT(N,D,E,M,W,z,IERR)                              
      IMPLICIT REAL*8(A-H,O-Z)                                          
c      INCLUDE 'parmar.f
      parameter(mdim=2**17,mr=30)
      INTEGER R,GROUP                                                   
      REAL*8  NORM,MACHEP                                               
      DIMENSION D(MDIM),E(MDIM),W(MDIM),RV1(MDIM),RV2(MDIM),   
     &            RV3(MDIM),RV4(MDIM),RV6(MDIM),z(mdim,mr)                         
c        COMMON A(MDIM,MDIM),Z(MDIM,MR)                                  
      MACHEP=2.0D0**(-37)                                               
      IERR=0                                                            
      DO 1 I=1,N-1                                                      
  1   E(I)=E(I+1)                                                       
      E(N)=0.0D0                                                        
      IF (M.EQ.0) GOTO 1001                                             
      ORDER=1.0D0                                                       
      DO 920 R=1,M                                                      
      ITS=1                                                             
      X1=W(R)                                                           
c	write(*,*)x1
      IF (R.NE.1) GOTO 100                                              
      NORM=0.0D0                                                        
      DO 40 I=1,N                                                       
  40  NORM=NORM+DABS(D(I))+DABS(E(I))                                   
c	write(*,*)norm
      EPS2=1.0D-3*NORM*DABS(ORDER)                                      
      write(*,*)eps2
      EPS3=MACHEP*NORM                                                  
c      write(*,*)eps3
      UK=DFLOAT(N)                                                      
      EPS4=UK*EPS3                                                      
      UK=EPS4/DSQRT(UK)                                                 
  80  GROUP=0                                                           
      GOTO 520                                                          
 100  IF (DABS(X1-X0).GE.EPS2) GOTO 80                                  
      GROUP=GROUP+1                                                     
      IF (ORDER*(X1-X0).LE.0.0D0) X1=X0+ORDER*EPS3                      
 520  V=0.0D0                                                           
      DO 580 I=1,N                                                      
      RV6(I)=UK                                                         
      IF (I.EQ.1) GOTO 560                                              
      IF (DABS(E(I)).LT.DABS(U)) GOTO 540                               
      XU=U/E(I-1)                                                       
      RV4(I)=XU                                                         
      RV1(I-1)=E(I-1)                                                   
      RV2(I-1)=D(I)-X1                                                  
      RV3(I-1)=0.0D0                                                    
      IF (I.NE.N) RV3(I-1)=E(I)                                         
      U=V-XU*RV2(I-1)                                                   
      V=-XU*RV3(I-1)                                                    
      GOTO 580                                                          
 540  XU=E(I-1)/U                                                       
c      write(*,*)e(i-1)
      RV4(I)=XU                                                         
      RV1(I-1)=U                                                        
      RV2(I-1)=V                                                        
      RV3(I-1)=0.0D0                                                    
 560  U=D(I)-X1-XU*V                                                    
      IF (I.NE.N) V=E(I)                                                
 580  CONTINUE                                                          
      IF (U.EQ.0.0) U=EPS3                                              
      RV1(N)=U                                                          
      RV2(N)=0.0D0                                                      
      RV3(N)=0.0D0                                                      
 600  DO 620 II=1,N                                                     
      I=1+N-II                                                          
      RV6(I)=(RV6(I)-U*RV2(I)-V*RV3(I))/RV1(I)                          
      V=U                                                               
 620  U=RV6(I)                                                          
      IF (GROUP.EQ.0) GOTO 700                                          
      DO 680 JJ=1,GROUP                                                 
      J=R-GROUP-1+JJ                                                    
      XU=0.0D0                                                          
      DO 640 I=1,N                                                      
 640  XU=XU+RV6(I)*Z(I,J)                                               
      DO 660 I=1,N                                                      
 660  RV6(I)=RV6(I)-XU*Z(I,J)                                           
 680  CONTINUE                                                          
 700  NORM=0.0D0                                                        
      DO 720 I=1,N                                                      
 720  NORM=NORM+DABS(RV6(I))                                            
      IF (NORM.GE.0.1D0) GOTO 840                                       
      IF (ITS.EQ.7) GOTO 830                                            
      ITS=ITS+1                                                         
      XU=EPS4/NORM                                                      
      DO 760 I=1,N                                                      
 760  RV6(I)=XU*RV6(I)                                                  
      DO 820 I=1,N                                                      
      U=RV6(I)                                                          
      IF (RV1(I-1).NE.E(I-1)) GOTO 800                                  
      U=RV6(I-1)                                                        
      RV6(I-1)=RV6(I)                                                   
 800  RV6(I)=U-RV4(I)*RV6(I-1)                                          
 820  CONTINUE                                                          
      GOTO 600                                                          
 830  IERR=-R                                                           
c       write(6,*)'  ITS=',ITS,'  R=',R,'  IERR=',IERR                  
      XU=0.0D0                                                          
      GOTO 870                                                          
 840  U=0.0D0                                                           
      DO 860 I=1,N                                                      
 860  U=U+RV6(I)**2                                                     
      XU=1.0D0/DSQRT(U)                                                 
 870  DO 900 I=1,N                                                      
      Z(I,R)=RV6(I)*XU
      IF (Z(1,R).LT.0.0D0) THEN
      Z(I,R)=-Z(I,R)
      ENDIF
 900  CONTINUE
      X0=X1                                                             
 920  CONTINUE                                                          
 1001 RETURN                                                            
      END                                                               
								
