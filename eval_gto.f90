subroutine eval_gto(ix,iy,iz,ao_value,ishell)
    use kinds, only: ikind,fp
    use reals, only: zero,one,three,onefive,seven,five
    use x2c1emem, only: nc,nunc
    use cvwint_module
    
    implicit none
    integer(kind=ikind) ishell,l,natm
    integer(kind=ikind) i,j,k
    integer(kind=ikind) nbas,nrc,istart,ipstart
    real(kind=fp) coord(3),dis(4)
    real(kind=fp) tmp,ix,iy,iz
    real(kind=fp),allocatable::prim(:),coeff(:),contr(:)
    real(kind=fp) ao_value(10*nc)
    
    nbas=int_nprims(ishell)
    nrc=int_nrc(ishell)
    istart=int_coff(ishell)
    ipstart=int_kprim(ishell)
    natm = int_katom(ishell)
    l = int_ktype(ishell)

    dis(1) = ix-int_coord(1,natm)
    dis(2) = iy-int_coord(2,natm)
    dis(3) = iz-int_coord(3,natm)
    dis(4) = dis(1)**2+dis(2)**2+dis(3)**2
  
    allocate(coeff(nbas*nrc),prim(nbas),contr(nrc))
    coeff = int_coef(istart+1:istart+nbas*nrc)
    prim = int_expnt(ipstart+1:ipstart+nbas) 
    
    do i = 1, nbas
        tmp = dexp(-prim(i)*dis(4))
        prim(i) = tmp
    end do     

    do j = 1, nrc
        tmp = zero
        do k = 1, nbas
            tmp = tmp +  coeff(k + (j - 1) * nbas) * prim(k)
        end do
            contr(j) = tmp      
    end do
 
    do j = 1, nrc
        select case(l)
        case(1)
            ao_value(j)=contr(j)
        case(2)
            ao_value(3*j-2) = contr(j) * dis(1)
            ao_value(3*j-1) = contr(j) * dis(2)
            ao_value(3*j) = contr(j) * dis(3)
        case(3)
            ao_value(6*j-5) = contr(j) * dis(1) * dis(1) * dsqrt(one/three)
            ao_value(6*j-4) = contr(j) * dis(1) * dis(2)
            ao_value(6*j-3) = contr(j) * dis(1) * dis(3) 
            ao_value(6*j-2) = contr(j) * dis(2) * dis(2) * dsqrt(one/three)
            ao_value(6*j-1) = contr(j) * dis(2) * dis(3)
            ao_value(6*j) = contr(j) * dis(3) * dis(3) * dsqrt(one/three)
        case(4)
            ao_value(10*j-9) = contr(j) * dis(1) * dis(1) * dis(1) * dsqrt(one/onefive)
            ao_value(10*j-8) = contr(j) * dis(1) * dis(1) * dis(2) * dsqrt(one/three)
            ao_value(10*j-7) = contr(j) * dis(1) * dis(1) * dis(3) * dsqrt(one/three)
            ao_value(10*j-6) = contr(j) * dis(1) * dis(2) * dis(2) * dsqrt(one/three)
            ao_value(10*j-5) = contr(j) * dis(1) * dis(2) * dis(3)
            ao_value(10*j-4) = contr(j) * dis(1) * dis(3) * dis(3) * dsqrt(one/three)
            ao_value(10*j-3) = contr(j) * dis(2) * dis(2) * dis(2) * dsqrt(one/onefive)
            ao_value(10*j-2) = contr(j) * dis(2) * dis(2) * dis(3) * dsqrt(one/three)
            ao_value(10*j-1) = contr(j) * dis(2) * dis(3) * dis(3) * dsqrt(one/three)
            ao_value(10*j) = contr(j) * dis(3) * dis(3) * dis(3) * dsqrt(one/onefive)
        case(5)
            contr(j) = contr(j) * dsqrt(one/three)
            ao_value(15*j-14) = contr(j) * dis(1) * dis(1) * dis(1) * dis(1) * dsqrt(one/seven) * dsqrt(one/five)  !xxxx
            ao_value(15*j-13) = contr(j) * dis(1) * dis(1) * dis(1) * dis(2) *  dsqrt(one/five)                    !xxxy
            ao_value(15*j-12) = contr(j) * dis(1) * dis(1) * dis(1) * dis(3) * dsqrt(one/five)                     !xxxz
            ao_value(15*j-11) = contr(j) * dis(1) * dis(1) * dis(2) * dis(2) * dsqrt(one/three)                    !xxyy
            ao_value(15*j-10) = contr(j) * dis(1) * dis(1) * dis(2) * dis(3)                                       !xxyz
            ao_value(15*j-9) = contr(j) * dis(1) * dis(1) * dis(3) * dis(3) * dsqrt(one/three)                     !xxzz
            ao_value(15*j-8) = contr(j) * dis(1) * dis(2) * dis(2) * dis(2) * dsqrt(one/five)                      !xyyy
            ao_value(15*j-7) = contr(j) * dis(1) * dis(2) * dis(2) * dis(3)                                        !xyyz 
            ao_value(15*j-6) = contr(j) * dis(1) * dis(2) * dis(3) * dis(3)                                        !xyzz
            ao_value(15*j-5) = contr(j) * dis(1) * dis(3) * dis(3) * dis(3) * dsqrt(one/five)                      !xzzz
            ao_value(15*j-4) = contr(j) * dis(2) * dis(2) * dis(2) * dis(2) * dsqrt(one/seven) * dsqrt(one/five)   !yyyy
            ao_value(15*j-3) = contr(j) * dis(2) * dis(2) * dis(2) * dis(3) * dsqrt(one/five)                      !yyyz
            ao_value(15*j-2) = contr(j) * dis(2) * dis(2) * dis(3) * dis(3)* dsqrt(one/three)                      !yyzz
            ao_value(15*j-1) = contr(j) * dis(2) * dis(3) * dis(3) * dis(3)* dsqrt(one/five)                       !yzzz
            ao_value(15*j) = contr(j) * dis(3) * dis(3) * dis(3) * dis(3) * dsqrt(one/seven) * dsqrt(one/five)     !zzzz
       case(6)
             ao_value(21*j-20) = contr(j) * dis(1) * dis(1) * dis(1) * dis(1) * dis(1) !xxxxx 
             ao_value(21*j-19) = contr(j) * dis(1) * dis(1) * dis(1) * dis(1) * dis(2) !xxxxy       
             ao_value(21*j-18) = contr(j) * dis(1) * dis(1) * dis(1) * dis(1) * dis(3) !xxxxz       
             ao_value(21*j-17) = contr(j) * dis(1) * dis(1) * dis(1) * dis(2) * dis(2) !xxxyy       
             ao_value(21*j-16) = contr(j) * dis(1) * dis(1) * dis(1) * dis(2) * dis(3) !xxxyz       
             ao_value(21*j-15) = contr(j) * dis(1) * dis(1) * dis(1) * dis(3) * dis(3) !xxxzz       
             ao_value(21*j-14) = contr(j) * dis(1) * dis(1) * dis(2) * dis(2) * dis(2) !xxyyy       
             ao_value(21*j-13) = contr(j) * dis(1) * dis(1) * dis(2) * dis(2) * dis(3) !xxyyz       
             ao_value(21*j-12) = contr(j) * dis(1) * dis(1) * dis(2) * dis(3) * dis(3) !xxyzz       
             ao_value(21*j-11) = contr(j) * dis(1) * dis(1) * dis(3) * dis(3) * dis(3) !xxzzz       
             ao_value(21*j-10) = contr(j) * dis(1) * dis(2) * dis(2) * dis(2) * dis(2) !xyyyy       
             ao_value(21*j-9) =  contr(j) * dis(1) * dis(2) * dis(2) * dis(2) * dis(3) !xyyyz       
             ao_value(21*j-8) =  contr(j) * dis(1) * dis(2) * dis(2) * dis(3) * dis(3) !xyyzz       
             ao_value(21*j-7) =  contr(j) * dis(1) * dis(2) * dis(3) * dis(3) * dis(3) !xyzzz       
             ao_value(21*j-6) =  contr(j) * dis(1) * dis(3) * dis(3) * dis(3) * dis(3) !xzzzz     
             ao_value(21*j-5) =  contr(j) * dis(2) * dis(2) * dis(2) * dis(2) * dis(2) !yyyyy        
             ao_value(21*j-4) =  contr(j) * dis(2) * dis(2) * dis(2) * dis(2) * dis(3) !yyyyz       
             ao_value(21*j-3) =  contr(j) * dis(2) * dis(2) * dis(2) * dis(3) * dis(3) !yyyzz       
             ao_value(21*j-2) =  contr(j) * dis(2) * dis(2) * dis(3) * dis(3) * dis(3) !yyzzz       
             ao_value(21*j-1) =  contr(j) * dis(2) * dis(3) * dis(3) * dis(3) * dis(3) !yzzzz       
             ao_value(21*j) =    contr(j) * dis(3) * dis(3) * dis(3) * dis(3) * dis(3) !zzzzz       
    end select
    end do    
    deallocate(coeff,prim,contr)
end subroutine eval_gto
