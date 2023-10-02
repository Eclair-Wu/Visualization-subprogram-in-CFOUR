subroutine evalao(core,maxcor)
      use kinds , only : ikind,fp
      implicit none
      integer(kind=ikind) :: maxcor
      real(kind=fp) :: core(maxcor)
!
      call intinfo
!
      call sox2c1eini(core,maxcor)
!
      call aogrid_drv(core,maxcor)

      call x2cint1ede 

end subroutine evalao
