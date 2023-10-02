subroutine findeig(eigen,a,d,aw,dw)
     
     use kinds,only:ikind,fp
     use x2c1emem,only:nmo,abovee

     implicit none 
     integer(kind=ikind) a,d,aw(nmo),dw(nmo)
     real(kind=fp) eigen(nmo)
     integer(kind=ikind) i

     
     a=0
     d=0
     
     do i = 1,nmo
        if (eigen(i) .gt. abovee) then
           a=a+1
           aw(a)=i
        else if (eigen(i) .lt. -abovee) then
           d=d+1
           dw(d)=i
        end if 
     end do        
      










end subroutine findeig
