subroutine genndo(eden)
      use kinds , only : ikind,fp,dcp
      use reals, only:zero,one
      use complexs , only :two_cp,zero_cp, one_cp
      use soccmod
      use soccmem
      use dg_module, only:dg_nbas
      implicit none
      integer(kind=ikind) i,j,l
      complex(kind=dcp) op1(2*dg_nbas),op2(3*dg_nbas-2)
      complex(kind=dcp) eden(dg_nbas,dg_nbas),gden(dg_nbas,dg_nbas),ndo(dg_nbas,dg_nbas)
      real(kind=fp) eig(dg_nbas)


     call rfile('symden',gden,dg_nbas*dg_nbas*2)
     do i = 1,no+ndrop
        gden(i,i)=one_cp+gden(i,i)
        eden(i,i)=one_cp+eden(i,i)
     end do

      ndo = eden-gden

      call zheev('V','U',dg_nbas,gden,dg_nbas,eig,op1,2*dg_nbas,op2,l)
      print*,'Ground state occ'
      call prvecr(eig,dg_nbas)
      call wfile('NOGtrans',gden,dg_nbas*dg_nbas*2)
      call wfile('NOGocc',eig,dg_nbas)
   
      call zheev('V','U',dg_nbas,eden,dg_nbas,eig,op1,2*dg_nbas,op2,l)
      print*,'Excited state occ'
      call prvecr(eig,dg_nbas)
      call wfile('NOEtrans',gden,dg_nbas*dg_nbas*2)
      call wfile('NOEocc',eig,dg_nbas)

      call zheev('V','U',dg_nbas,ndo,dg_nbas,eig,op1,2*dg_nbas,op2,l)
      print*,'NDO eigenvalue'
      call prvecr(eig,dg_nbas)
      call wfile('NDOocc',eig,dg_nbas)
      call wfile('NDOtrans',ndo,dg_nbas*dg_nbas*2)

     
end subroutine genndo
