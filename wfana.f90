subroutine wfana(work)
     use x2c1emem
     use kinds,only:ikind,dcp,fp
     use complexs,only:one_cp,zero_cp

     real(kind=fp) work(lwork)
     real(kind=fp) eig(nmo)
     complex(kind=dcp) mo_occ(2*nc,nmo),st(2*nc),no_trans(nmo,nmo),no(2*nc,nmo)
     character(len=20) vtype,cubet 
     character(len=5) order
     integer(kind=ikind) i,j,l,anum,dnum,awhere(nmo),dwhere(nmo)
     complex(kind=dcp) htrans(occ,occ),ptrans(unocc,unocc),ntoh(2*nc,occ),ntop(2*nc,unocc)
     real(kind=fp),allocatable::sig(:) 


     call rfile('scf2ccmo',mo_occ,nmo*4*nc)


     select case(anatype)

     case(0)
          st = mo_occ(:,state)
          print*,"Visualization a given HF MO"
          select case(state)
            case(1:9)
              write(order,'(I1)') state
            case(10:99)
              write(order,'(I2)') state
          end select
          vtype='HF'
          write(cubet,*) trim(vtype),order
          call eval_mo(work(aovalue),st,cubet)   
     case (1)     
         print*,"Visualization NDO after EOMCCSD"
         call rfile('NDOocc',eig,nmo)
         call rfile('NDOtrans',no_trans,nmo*nmo*2)
         call zgemm('N','N',2*nc,nmo,2*nc,one_cp,mo_occ,2*nc,no_trans,nmo,zero_cp,no,2*nc)
         call findeig(eig,anum,dnum,awhere,dwhere)
         do i = 1,dnum
            j = dwhere(i)
            st=no(:,j)
            write(order,'(I1)') i
            vtype='NDO_D'
            write(cubet,*) trim(vtype),order
            print*,'Detachment occ:',eig(j),cubet
         call eval_mo(work(aovalue),st,cubet)
         end do
         do i = 1,anum
            j = awhere(i)
            st=no(:,j)
            write(order,'(I1)') i
            vtype='NDO_A'
            write(cubet,*) trim(vtype),order
            print*,'Attachment occ:',eig(j),cubet

            call eval_mo(work(aovalue),st,cubet)
         end do
     case(2)
         print*,"Visualization NO after CCSD"
         call rfile('NOGtrans',no_trans,nmo*nmo*2)
         call rfile('NOGocc',eig,nmo)
         call zgemm('N','N',2*nc,nmo,2*nc,one_cp,mo_occ,2*nc,no_trans,nmo,zero_cp,no,2*nc)
         call findeig(eig,anum,dnum,awhere,dwhere)
          select case(state)
            case(1:9)
              write(order,'(I1)') state
            case(10:99)
              write(order,'(I2)') state
          end select      
         vtype='NO_G'
         st = no(:,nmo+1-state)
         print*,'Occ of this NO:',eig(nmo+1-state)
         write(cubet,*) trim(vtype),order
         call eval_mo(work(aovalue),st,cubet)
     case(3)
         print*,"Visualization NO after EOMCCSD"
         call rfile('NOEtrans',no_trans,nmo*nmo*2)
         call rfile('NOEocc',eig,nmo)
         call zgemm('N','N',2*nc,nmo,2*nc,one_cp,mo_occ,2*nc,no_trans,nmo,zero_cp,no,2*nc)
         call findeig(eig,anum,dnum,awhere,dwhere)
         select case(state)
            case(1:9)
              write(order,'(I1)') state
            case(10:99)
              write(order,'(I2)') state
          end select
         st = no(:,nmo+1-state)
         print*,'Occ of this NO:',eig(nmo+1-state)
         vtype='NO_E'
         write(cubet,*) trim(vtype),order
         call eval_mo(work(aovalue),st,cubet)
    case(4)
         if (occ.ge.unocc) then
            l = unocc
         else
            l = occ
         end if
         allocate(sig(l))
         print*,"Visualization NTO after EOMCCSD"  
         call rfile('NTO_h',htrans,occ*occ*2)
         call rfile('NTOocc',sig,l)
         call rfile('NTO_p',ptrans,unocc*unocc*2)
         call zgemm('N','C',2*nc,occ,occ,one_cp,mo_occ(:,ndrop+1:ndrop+occ),2*nc,htrans,occ,zero_cp,ntoh,2*nc)
         call zgemm('N','N',2*nc,unocc,unocc,one_cp,mo_occ(:,occ+ndrop+1:),2*nc,ptrans,unocc,zero_cp,ntop,2*nc)
         call findeig(sig,anum,dnum,awhere,dwhere)
         call prvecr(sig,l)
         print*,occ,unocc,ndrop
         do i = 1,anum
            j = awhere(i)
            st=ntoh(:,j)
            write(order,'(I1)') i
            vtype='NTO_H'
            write(cubet,*) trim(vtype),order
            print*,'Weights of this NTO pair:',sig(j)**2,i
            call eval_mo(work(aovalue),st,cubet)
            st=ntop(:,j)
            vtype='NTO_P'
            write(cubet,*) trim(vtype),order
            call eval_mo(work(aovalue),st,cubet)
         end do
         deallocate(sig)
     end select





end subroutine wfana
