
 module mod_calc_vol_cv_centers
  contains
      subroutine calc_vol_cv_centers_uns(xc,yc,zc,vol,x,y,z,e2vx,ef2nb,ef2nb_idx,etype,nvx,nelem,nbndry,nf,ne2vxmax,nsec,esec,xyzip,aip)
      use mod_util
      implicit none

      integer :: nvx,nelem,nf,ne2vxmax,nbndry,nsec,esec(2,nsec)
      integer :: e2vx(ne2vxmax,nelem),etype(nelem),ef2nb(2*nf-nbndry,2),ef2nb_idx(nelem-nbndry+1)
      real :: xc(nelem),yc(nelem),zc(nelem),x(nvx),y(nvx),z(nvx),xyzip(3,nf),aip(3,nf),vol(nelem-nbndry)
!
      integer :: e,s,v,t,nv,idx,gf,lf,enb,gf_sgn
      real :: cs_cntr(3),geom_cntr(3),h(3),sub_vol,sub_cntr(3),sum_vol,cs_area,rc(3)
!
      do s=1,nsec
        do e=esec(1,s),esec(2,s)
          if(e<=nelem-nbndry) then ! 3d
            t=etype(s)

            nv=element_nvx(t)
            geom_cntr = 0.
            do idx=1,nv
              v=e2vx(idx,e)
              geom_cntr = geom_cntr + [x(v),y(v),z(v)]
            enddo
            geom_cntr=geom_cntr/nv

            sum_vol=0.
            sub_cntr=0.
            rc=0.
            do idx=ef2nb_idx(e),ef2nb_idx(e+1)-1
              gf=ef2nb(idx,2)
              gf_sgn=sgn(gf)
              gf=abs(gf)
              cs_cntr = xyzip(:,gf)
              !!! set the center of halo if neighbor is boundary
              call get_idx(ef2nb(idx,1),0,enb,lf)
              if(lf==0) then
                xc(enb)=xyzip(1,gf)
                yc(enb)=xyzip(2,gf)
                zc(enb)=xyzip(3,gf)
              end if
              !!!
              cs_area = sqrt(aip(1,gf)**2+aip(2,gf)**2+aip(3,gf)**2)
              h=cs_cntr-geom_cntr
              sub_vol=dot_product(gf_sgn*aip(:,gf),h)/(3.)
              sub_cntr= (0.25)*geom_cntr+(0.75)*cs_cntr
              sum_vol=sum_vol+sub_vol
              rc=rc+sub_cntr*sub_vol
            end do
            rc=rc/sum_vol
            xc(e)=rc(1)
            yc(e)=rc(2)
            zc(e)=rc(3)
            vol(e)=sum_vol
            if(vol(e) <= 1e-20) write(*,*) 'Negative or zero volume detected @ Location',rc,'Volume',vol(e)
          endif
        enddo
      enddo

      end subroutine

  end module

