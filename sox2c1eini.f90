subroutine sox2c1eini(work,length)  
! 
      use kinds , only : ikind,fp,dcp
      use reals , only : one,zero,two,three,five
      use x2c1emem
      use x2cint1emd
      use cvwint_module
!
      implicit none
!
      real(kind=fp) :: work(length),x,y,z,r
      integer(kind=ikind) n,lworktmp,length,i,j,cptran,eof,tmp
      character*8 :: string
!      
      call x2cint1eini
!      
      lwork=length
      n=nunc
      ndimsq=n*n*16

!
      cptran=1
      somf=cptran+10000
      iend=somf+n*n*32
!
      
      ill=iend
      ils=ill+n*n*8
      isl=ils+n*n*8
      iss=isl+n*n*8
      iend=iss+n*n*8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      t=iend
      v=t  +ndim1e
      s=v  +ndim1e
      pvp=s  +ndim1e
      paopvpso=pvp+ndim1e
      iend=paopvpso+3*ndim1e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      c=iend
      cl=c+2*ndimsq
      cs=cl+n*n*8
      cml=cs+n*n*8
      cms=cml+n*n*8
      iscr1=cms+n*n*8
      iscr2 =iscr1 +2*16*n*n
      iscr3 =iscr2 +2*16*n*n
      iscr4 =iscr3 +2*16*n*n
      iscr5 =iscr4 +2*16*n*n
      iscr6 =iscr5 +2*16*n*n
      iscr7 =iscr6 +2*16*n*n
      iscr8 =iscr7 +2*16*n*n
      iscr9 =iscr8 +2*16*n*n
      iscr10=iscr9 +2*16*n*n
      iscr11=iscr10+2*16*n*n
      iscr12=iscr11+2*16*n*n
      iend  =iscr12+2*ndimsq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      h4c =iend
      c4c =h4c+32*n*n
      s4c =c4c+32*n*n
      iend=s4c+32*n*n
!
      ifac   = iend
      ifac2  = ifac   +    n*8
      ifacl  = ifac2  +    n*8
      ifaclc = ifacl  +    n*8
      ieig   = ifaclc +    n*8
      ihfw   = ieig   +    n*8
      ix     = ihfw   +  n*n*8
      ir     = ix     +  n*n*8
      il     = ir     +  n*n*8
      ist    = il     +  n*n*8
      ssq    = ist    +  n*n*8
      ssqi   = ssq    +  n*n*8
      iend   = ssqi   +  n*n*8
!
      seig   = iend
      svec   = seig   +    n*4
      amat   = svec   +  n*n*8
      aeig   = amat   +  n*n*8
      avec   = aeig   +    n*4
      asq    = avec   +  n*n*8
      asqi   = asq    +  n*n*8
      xatmp  = asqi   +  n*n*8
      iend   = xatmp  +  n*n*8
!
      v1     = iend
      v1a    = v1     +  n*n
      v1b    = v1a    +  n*n
      v1c    = v1b    +  n*n
      ovlps  = v1c    +  n*n
      ovlpss = ovlps  +  n*n
      ssq1c  = ovlpss +  n*n
      ssqi1c = ssq1c  +  n*n
      pvp1   = ssqi1c +  n*n
      sopvp1 = pvp1   +  n*n
      s1ori  = sopvp1 +  n*n*3
      t1     = s1ori  +  n*n
      ihll1  = t1     +  n*n  
      ihls1  = ihll1  +  n*n*8
      ihsl1  = ihls1  +  n*n*8
      ihss1  = ihsl1  +  n*n*8
      isll0  = ihss1  +  n*n*8
      isss0  = isll0  +  n*n*8
      ssqchi = isss0  +  n*n*8
      ssqichi= ssqchi +  n*n*8
      isll1  = ssqichi+  n*n*8
      isss1  = isll1  +  n*n*8
      v14c   = isss1  +  n*n*8
      s14c   = v14c   +  n*n*32
      xx     = s14c   +  n*n*32
      xy     = xx     +  n*n
      xz     = xy     +  n*n
      yx     = xz     +  n*n
      yy     = yx     +  n*n
      yz     = yy     +  n*n
      zx     = yz     +  n*n
      zy     = zx     +  n*n
      zz     = zy     +  n*n
      zxx    = zz     +  n*n
      xzx    = zxx    +  n*n
      yzx    = xzx    +  n*n
      xyz    = yzx    +  n*n
      yyz    = xyz    +  n*n
      zyy    = yyz    +  n*n
      zzx    = zyy    +  n*n
      xzz    = zzx    +  n*n
      yzz    = xzz    +  n*n
      zyz    = yzz    +  n*n
      xzy    = zyz    +  n*n
      yzy    = xzy    +  n*n
      iend   = yzy    +  n*n
!      
      geogrd = iend
      iend   =geogrd  +  n*n*8*int_natms*6
!      
      c1     = iend 
      cl1    = c1     +  n*n*32
      cs1    = cl1    +  n*n*8
      ihfw1  = cs1    +  n*n*8
      l1     = ihfw1  +  n*n*8
      x1     = l1     +  n*n*8
      r1     = x1     +  n*n*8
      iend   = r1     +  n*n*8
!
      grids=0
      interval=0.0d0
      tmp=0
      eof=0
      open(unit=30,file='ZMAT',form='formatted')
      rewind(30)
      do while(tmp.eq.0 .and. eof.eq.0)
       read(30,'(A)',iostat=eof) string
       tmp = index(string,"%aogrid")
       if(tmp.ne.0) read(30,*)grids,interval,plottype
      end do
      
      tmp=0
      eof=0
      close(unit=30,status='keep')
      open(unit=30,file='ZMAT',form='formatted')
      rewind(30)
      do while(tmp.eq.0 .and. eof.eq.0)
       read(30,'(A)',iostat=eof) string
       tmp = index(string,"%visual")
       if(tmp.ne.0) read(30,*)anatype,state,abovee
      end do
      close(unit=30,status='keep')
     
      
      ng=2*grids+1

      aovalue=iend
      iend=aovalue+ng*ng*ng*nc
      call rfile('nmototal',nmo,1) 
      call rfile('unoccmo',unocc,1)
      call rfile('occmo',occ,1)
      
      ndrop=nmo-unocc-occ
     
      call get_trafo(work(cptran),work(iscr1),3)

end subroutine sox2c1eini
