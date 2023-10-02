subroutine eval_mo(ao_value,st,vtype)
      use kinds, only:ikind,fp
      use x2c1emem, only:nmo,ng,nc,state,plottype
      use reals, only:one,zero

      integer(kind=ikind) i,j, ngrid 
      real(kind=fp) st(4*nc)  
      real(kind=fp),allocatable::mo_den(:),mo_reala(:),mo_imaga(:),mo_r(:),mo_realb(:),mo_imagb(:)
      real(kind=fp)::ao_value(ng**3*nc)
      real(kind=fp),allocatable::mo_dena(:),mo_denb(:),mo_a(:),mo_b(:),mo_imag(:),mo_real(:)
      character(len=20) visual(8),vtype

      visual(1) = 'den.cube'
      visual(2) = 'den_a.cube'
      visual(3) = 'den_b.cube'
      visual(4) = 'mo.cube'
      visual(5) = 'mo_a.cube'
      visual(6) = 'mo_b.cube' 
      visual(7) = 'mo_real.cube'
      visual(8) = 'mo_imag.cube'

      ngrid=ng**3

      
     
      allocate(mo_realb(ngrid),mo_imagb(ngrid),mo_reala(ngrid),mo_imaga(ngrid),mo_r(ngrid),mo_den(ngrid))
      allocate(mo_dena(ngrid),mo_denb(ngrid),mo_a(ngrid),mo_b(ngrid),mo_real(ngrid),mo_imag(ngrid))
      do i = 1, ngrid
          mo_reala(i) = zero
          mo_imaga(i) = zero
          mo_realb(i) = zero
          mo_imagb(i) = zero
         do j = 1,nc
          mo_reala(i) = mo_reala(i)+ao_value((i-1)*nc+j)*st(2*j-1)
          mo_imaga(i) = mo_imaga(i)+ao_value((i-1)*nc+j)*st(2*j)
          mo_realb(i) = mo_realb(i)+ao_value((i-1)*nc+j)*st(2*nc+2*j-1)
          mo_imagb(i) = mo_imagb(i)+ao_value((i-1)*nc+j)*st(2*nc+2*j)
         end do
          mo_dena(i)=mo_reala(i)**2+mo_imaga(i)**2
          mo_denb(i)=mo_realb(i)**2+mo_imagb(i)**2
          mo_den(i)=mo_dena(i)+mo_denb(i)
          mo_r(i)=dsqrt(mo_den(i))
      end do

      select case(plottype)
      case(0)
          print*,'Generate density and module of mo'
          call gencube(mo_den,visual(1),vtype)
          call gencube(mo_r,visual(4),vtype)
      case(1)
         print*,'Generate alpha and beta density '     
         call gencube(mo_dena,visual(2),vtype)
         call gencube(mo_denb,visual(3),vtype)
      case(2)
         print*,'Generate module of alpha and beta mo' 
         do i = 1,ngrid
            mo_a(i)=dsqrt(mo_dena(i))
            mo_b(i)=dsqrt(mo_denb(i))
         end do   

         call gencube(mo_a,visual(5),vtype)
         call gencube(mo_b,visual(6),vtype)
      case(3)
         print*,'Generate alpha and beta density,module of alpha and beta mo' 
         do i = 1,ngrid
            mo_a(i)=dsqrt(mo_dena(i))
            mo_b(i)=dsqrt(mo_denb(i))
         end do
         call gencube(mo_dena,visual(5),vtype)
         call gencube(mo_denb,visual(6),vtype)
         call gencube(mo_a,visual(2),vtype)
         call gencube(mo_b,visual(3),vtype)
      case(4)
         print*,'Generate real and imaginary part of mo' 
         do i = 1,ngrid
             mo_real(i)=mo_reala(i)+mo_realb(i)
             mo_imag(i)=mo_imaga(i)+mo_imagb(i)
         end do
         call gencube(mo_real,visual(7),vtype)
         call gencube(mo_imag,visual(8),vtype)

      end select






      
      deallocate (mo_real,mo_imag,mo_reala,mo_imaga,mo_den,mo_realb,mo_imagb,mo_r,mo_dena,mo_denb,mo_a,mo_b)
      end subroutine eval_mo
