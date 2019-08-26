        implicit real*8(a-h,o-z)
        parameter (ntp=10000)
        dimension vlevel(ntp),xo(ntp),xn(ntp)
        dimension vl(ntp),xlo(ntp),xln(ntp)
        dimension ni(ntp),nj(ntp),xxi(ntp),xxj(ntp)
        open(1,file='mgh2-cavnz-cvnz-ciq-j0.dat')
        open(2,file='podvr-2.out')
        read(2,*)
        read(2,*)
          k=0
        do 100 i=1,ntp
        read(2,*)ni(i),vl(i),xlo(i),xln(i),xxi(i)
        if((abs(vl(i)-vl(i-1)).gt.0.001).and.
     &   (dabs(xlo(i)).gt.0.00001).and.(dabs(xln(i)).lt.0.000001))then
         write(1,999)ni(i),vl(i),xlo(i),xln(i),xxi(i)
         endif 

        if(ni(i).gt.6000) then
        goto 300 
        endif 

 
        do 200 j=1,200
        read(2,*)nj(j),vlevel(j),xo(j),xn(j),xxj(j)
        if((abs(vlevel(j)-vl(i)).gt.0.001).and.
     &  (dabs(xo(j)).gt.0.00001).and.(dabs(xn(j)).lt.0.000001))then
        write(1,999)nj(j),vlevel(j),xo(j),xn(j),xxj(j)
        go to 100
        endif 
        if((dabs(vlevel(j)-vl(i)).lt.0.0001).and.
     &  (dabs(xo(j)).gt.0.00001).and.(dabs(xn(j)).lt.0.000001))then 
        write(1,999)nj(j),vlevel(j),xo(j),xn(j),xxj(j)
        endif
200     continue
100     continue       
300     stop 
999     format(1x,i7,f20.12,10e15.6)
        end



