subroutine gencube(array,cubetype,vtype)

      use x2c1emem
      use kinds, only:ikind,fp
      use cvwint_module
      use reals, only:zero,three,one
      implicit none
      integer(kind=ikind) i,j,k,l
      real(kind=fp)::array(ng**3)
      real(kind=fp) box    
      character(len=20) cubetype,vtype
      character(len=40) title

      write(title,*) trim(vtype),trim(cubetype)
      
      open(1,file = title,status = 'replace')
      write(1,*) 'MO/Electron density in real space (e/Bohr^3)'
      write(1,*) 'Written by Yixuan, just a test'
      box = interval*grids
      write(1,100) int_natms,-box,-box,-box
      write(1,100) ng,interval,zero,zero
      write(1,100) ng,zero,interval,zero
      write(1,100) ng,zero,zero,interval
      do i = 1,int_natms
        write(1,200) Int(int_charge(i)),zero,int_coord(1,i),int_coord(2,i),int_coord(3,i)
      end do
      do i  = 1,ng
        do j = 1,ng
          do k  = 1,ng
            if (MOD(k,6).eq.0) then
              write(1,400) array((i-1)*ng*ng+(j-1)*ng+k)
            else if(k.eq.ng) then
              write(1,400) array((i-1)*ng*ng+(j-1)*ng+k)
            else 
              write(1,300) array((i-1)*ng*ng+(j-1)*ng+k)
            end if
          end do
        end do
      end do
      close(1)
      


100 format(I5,F12.6,F12.6,F12.6)
200 format(I5,F12.6,F12.6,F12.6,F12.6)
300 format(F12.6,' ',\)
400 format(F12.6)
      end subroutine gencube     
