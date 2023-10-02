subroutine intinfo
      use cvwint_module
      use basinfo2_module
      use basinfo1_module
      use x2c1emem
      implicit none
!
      call int_itrf0
      call int_itrf1
      call int_itrf3
      call int_itrffg
      call int_itrfg
      call int_itrfh
!
      call int_itrf4
      nunc=int_nsph
      nc=int1_nsph
      ndim1e=nunc*(nunc+1)/2
      call int_itrf2
!
end subroutine intinfo

