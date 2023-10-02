subroutine gennto(nto)
      use kinds , only : ikind,fp,dcp
      use reals, only:zero
      use complexs , only :two_cp,zero_cp, one_cp
      use soccmod
      use soccmem
      use dg_module, only:dg_nbas
      implicit none
      integer(kind=ikind) i,j,l
      complex(kind=dcp) op1(nu,nu),op2(no,no),op3(3*dg_nbas),op4(5*dg_nbas)
      complex(kind=dcp) nto(nu,no)
      real(kind=fp),allocatable::sig(:)
      call wfile('occmo',no,1)
      call wfile('unoccmo',nu,1)
      
      if (nu.ge.no) then
         l=no
      else
         l=nu
      end if
      allocate(sig(l))
     
      call zgesvd('A','A',nu,no,nto,nu,sig,op1,nu,op2,no,op3,3*dg_nbas,op4,i)
      call wfile('NTOocc',sig,l)     
      call wfile('NTO_p',op1,nu*nu*2)
      call wfile('NTO_h',op2,no*no*2)
      print*,'NTO eigenvalue'
      call prvecr(sig,l)
      deallocate (sig)
end subroutine gennto
