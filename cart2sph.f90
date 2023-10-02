subroutine cart2sph(ao_value,ishell,ao_sph)
      use kinds,only:ikind,fp
      use reals
      use x2c1emem,only:nc,nunc
      use cvwint_module
      implicit none

      integer(kind=ikind) ishell,l,i,nrc
      real(kind=fp) ao_value(10*nc),ao_sph(7*nc)

      l=int_ktype(ishell)
      nrc=int_nrc(ishell)
      select case(l)
      case(1)
         do i = 1,nrc
         ao_sph(i) = ao_value(i)
         end do
      case(2)
         do i = 1,nrc
         ao_sph(i)=ao_value(i*3-2)
         ao_sph(nrc+i)=ao_value(i*3-1)
         ao_sph(nrc*2+i)=ao_value(i*3)
         end do 
      case(3)
         do i = 1,nrc
            ao_sph(4*nrc+i) = ao_value(6*i)-(ao_value(6*i-5)+ao_value(6*i-2))/two  !lz=0 
            ao_sph(nrc+i) = ao_value(6*i-4)  !lz=-2
            ao_sph(2*nrc+i) = ao_value(6*i-3)  !lz=+1
            ao_sph(i) = (ao_value(6*i-5)-ao_value(6*i-2))*dsqrt(three/four)  !lz=+2
            ao_sph(3*nrc+i) = ao_value(6*i-1)  !lz=-1
         end do
      case(4)
         do i = 1,nrc
            ao_sph(i+4*nrc) = ao_value(10*i-4)*dsqrt(six/five)-ao_value(10*i-9)*dsqrt(six)/four-ao_value(10*i-6)*dsqrt(six/five)/four  !lz=+1
            ao_sph(i+5*nrc) = ao_value(10*i-1)*dsqrt(six/five)-ao_value(10*i-3)*dsqrt(six)/four-ao_value(10*i-8)*dsqrt(six/five)/four  !lz=-1
            ao_sph(i+6*nrc) = ao_value(10*i)-(ao_value(10*i-7)+ao_value(10*i-2))*dsqrt(nine/twenty)  !lz=0
            ao_sph(i) = ao_value(10*i-9)*dsqrt(five/eight)-ao_value(10*i-6)*dsqrt(nine/eight)  !lz=+3
            ao_sph(i+3*nrc) = ao_value(10*i-5)  !lz=-2
            ao_sph(i+nrc) = ao_value(10*i-8)*dsqrt(nine/eight)-ao_value(10*i-3)*dsqrt(five/eight)  !lz=-3
            ao_sph(i+2*nrc) = (ao_value(10*i-7)-ao_value(10*i-2))*dsqrt(three/four)  !lz=+2
        end do
      case(5)
         do i = 1,nrc
            ao_sph(nrc*8+i) =ao_value(15*i)+(ao_value(15*i-4)+ao_value(15*i-14))*three/eight-(ao_value(15*i-9)+ao_value(15*i-2)-ao_value(15*i-11)/four)*dsqrt(nine/seven)*dsqrt(three/five)  !lz=0
            ao_sph(nrc*5+i) = ao_value(15*i-6)*dsqrt(nine/seven)-(ao_value(15*i-13)+ao_value(15*i-8))*dsqrt(five/seven)/two  !lz=-2
            ao_sph(nrc*6+i) = (ao_value(15*i-5)-(three/four)*ao_value(15*i-12))*dsqrt(ten/seven)-(ao_value(15*i-7))*dsqrt(two/seven)*(three/four)  !lz=+1
            ao_sph(i) = (ao_value(15*i-4)+ao_value(15*i-14))*dsqrt(five*seven)/eight-ao_value(15*i-11)*dsqrt(three)*three/four  !lz=+4
            ao_sph(nrc*3+i) = ao_value(15*i-10)*dsqrt(two)*three/four-ao_value(15*i-3)*dsqrt(ten)/four  !lz=-3
            ao_sph(nrc*4+i) = (ao_value(15*i-9)-ao_value(15*i-2))*dsqrt(three/seven)*three/two-(ao_value(15*i-14)-ao_value(15*i-4))*dsqrt(five)/four  !lz=+2
            ao_sph(nrc+i) = (ao_value(15*i-13)-ao_value(15*i-8))*dsqrt(five)/two  !lz=-4
            ao_sph(nrc*2+i) = ao_value(15*i-12)*dsqrt(ten)/four-ao_value(15*i-7)*dsqrt(two)*three/four  !lz=+3
            ao_sph(nrc*7+i) = ao_value(15*i-6)*three/dsqrt(seven)-(ao_value(15*i-13)-ao_value(15*i-8))*dsqrt(five/seven)/two   !lz=-1
        end do
     case(6)
         do i = 1,nrc
         ao_sph(i)=zero
         ao_sph(nrc+i)=zero
         ao_sph(nrc*2+i)=zero
         ao_sph(nrc*3+i)=zero      
         ao_sph(nrc*4+i)=zero
         ao_sph(nrc*5+i)=zero
         ao_sph(nrc*6+i)=zero
         ao_sph(nrc*7+i)=zero
         ao_sph(nrc*8+i)=zero
         ao_sph(nrc*9+i)=zero
         ao_sph(nrc*10+i)=zero
        end do
     end select
end subroutine cart2sph


