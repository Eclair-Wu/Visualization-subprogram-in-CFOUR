subroutine aogrid_drv(work,lwork)
      use kinds, only: ikind,fp
      use reals, only: one
      use cvwint_module
      use x2c1emem,only:grids,aovalue,interval,ng,nc
      implicit none

      integer(kind=ikind) ishell,lshell,i,j,k,lwork
      integer(kind=ikind) nrc,ndeg,iidim,ndim,iskip,ndegsph,ndegcart
      real(kind=fp) work(lwork)
      real(kind=fp) x,y,z
      real(kind=fp),allocatable:: ao_value(:),aoval(:),aosph(:)
      logical yesno

      iskip=ng**3
      ndim=1
     inquire(file='AO_GRID',exist=yesno)
     if (yesno) then
        call rfile('AO_GRID',work(aovalue),ng**3*nc)
        call wfana(work)
     else  
     allocate(ao_value(ng**3*nc))
      do i=-grids,grids
       x=i*interval
       do j=-grids,grids
        y=j*interval
        do k=-grids,grids
         z=k*interval
         do ishell=1,int_nshell
          lshell=int_ktype(ishell)
          ndegsph=2*lshell-1
          nrc=int_nrc(ishell)
          ndegcart=lshell*(lshell+1)/2
          allocate(aoval(nrc*ndegcart),aosph(nrc*ndegsph))
          call eval_gto(x,y,z,aoval,ishell)
          call cart2sph(aoval,ishell,aosph)
          call dcopy(nrc*ndegsph,aosph,1,ao_value(ndim:ndim+ndegsph),1)
          ndim = ndim+ndegsph*nrc
          deallocate(aoval,aosph)
         end do
        end do
       end do
      end do
      call dcopy(ng**3*nc,ao_value,1,work(aovalue),1)
      call wfile('AO_GRID',work(aovalue),ng**3*nc)
      call wfana(work)
      deallocate(ao_value)
     end if
end subroutine aogrid_drv
