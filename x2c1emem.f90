module x2c1emem  
      use kinds , only : ikind,fp
      implicit none
!
      integer(kind=ikind) :: t,v,s,pvp,paopvpso,ieig,somf
!
      integer(kind=ikind) :: ndimsq,ndim1e
!
      integer(kind=ikind) :: c,cl,cs,cml,cms     
!
      integer(kind=ikind) :: iscr1,iscr2,iscr3,iscr4,iscr5,iscr6
      integer(kind=ikind) :: iscr7,iscr8,iscr9,iscr10,iscr11,iscr12
!
      integer(kind=ikind) :: lwork,iend,sleft,istart,iendbp
!
      integer(kind=ikind) :: h4c,c4c,s4c,ifac,ifac2,ihfw,ix,ir,ist
      integer(kind=ikind) :: ill,ils,isl,iss,il
!
      integer(kind=ikind) :: v1,pvp1,sopvp1,x1,r1,l1,ihfw1
      integer(kind=ikind) :: v1a,v1b,v1c
      integer(kind=ikind) :: cl0,sss,sst
      integer(kind=ikind) :: ssq,ssqi
      integer(kind=ikind) :: seig,svec,amat,aeig,avec,asq,asqi,xatmp
      integer(kind=ikind) :: ihll0,ihls0,ihsl0,ihss0
      integer(kind=ikind) :: ihll1,ihls1,ihsl1,ihss1
      integer(kind=ikind) :: isll0,isll1,isss0,isss1
      integer(kind=ikind) :: ssqchi,ssqichi
      integer(kind=ikind) :: v14c ,s14c, cl1  ,cs1  ,c1, s1ori, t1
      integer(kind=ikind) :: ovlps,ovlpss,ssq1c,ssqi1c
!
      integer(kind=ikind) :: ix1flg,nc,nunc,iatmx
      integer(kind=ikind) :: xx,xy,xz,yx,yy,yz,zx,zy,zz      
      integer(kind=ikind) :: zxx,xzx,yzx,xyz,yyz,zyy,zzx,xzz,yzz,zyz,xzy,yzy
      integer(kind=ikind) :: ifacl,ifaclc
!      
      integer(kind=ikind) :: geogrd,occ,unocc,ndrop
      integer(kind=ikind) :: grids,ng,state,aovalue,nmo,anatype,plottype
      !ao_value,aoval,ao_sph
      real(kind=fp)       :: interval,abovee
end module x2c1emem
