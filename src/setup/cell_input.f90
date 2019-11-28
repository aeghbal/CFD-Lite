!
!  cell_input.f90
!
    module mod_cell
      use mod_mg_lvl
      implicit none

    contains

    subroutine cell_input(geom,cgfilename)
    use mod_calc_vol_cv_centers
    use mod_calc_aip_xyzip
    use mod_agglomeration
    use mod_cgns
    use mod_util
    implicit none
    character (len=32) :: agg_method
    character (len=180) :: cgfilename
    character (len=64) :: cval

    type (cgnsdb) :: cg
    type (geometry_t) :: geom

    integer :: ne,nf,nbf,nsec,target_ncv,aiplength,length,vlength,p,istat
    integer :: m,n,g,gnb,ls,idx,s,is,ie,js,je,ks,ke,i,j,k,np,cvs(3),cvs_f(3),l,step(3),z,gr(3),gf2vx_size
    integer , pointer , dimension(:) :: nl,ne2vx
    integer ::icell,parent_flag,total_elemcnt,total_datasize,sec
    integer :: nvx,nelem,nbndry,ne2vxmax,npmax,range(6),mem_size(2),vx2e_size,face
    real, parameter :: unitscale=1.0_8

     icell= 1
     agg_method = 'Weighted Directional'!'Oct-Tree Isotropic'

     associate(mg_lvl => geom%mg)

       cg%file=cgfilename
       call cgns_db_open(cg)
       call cg_zone_read_f(cg%unit, cg%base, 1, mg_lvl%cellname, cvs, ier)
       nvx=cvs(1)

       call cg_goto_f(cg%unit, cg%base, ier, 'Zone_t', icell, 'Elements_t', 1, 'end')
       if (ier .eq. ERROR) call cg_error_exit_f
       call cg_nsections_f(cg%unit, cg%base, icell, nsec, ier)
       if (ier .eq. ERROR) call cg_error_exit_f


       allocate(mg_lvl%esec(2,nsec))
       allocate(mg_lvl%etype(nsec))
       allocate(mg_lvl%ne2vx(nsec))
       allocate(mg_lvl%sectionName(nsec))

       cvs = 0
       ne2vxmax = 0
       vx2e_size = 0
       do sec=1,nsec
         call cg_section_read_f(cg%unit, cg%base, icell,sec,mg_lvl%sectionName(sec),mg_lvl%etype(sec),&
                                    mg_lvl%esec(1,sec),mg_lvl%esec(2,sec),nbndry,parent_flag, ier)
         if (ier .eq. ERROR) call cg_error_exit_f
         if(mg_lvl%etype(sec)>=10 .and. mg_lvl%etype(sec)<=20) then ! 3d cells
           cvs(1) = cvs(1) + mg_lvl%esec(2,sec)-mg_lvl%esec(1,sec)+1
           cvs(2) = cvs(2) + element_nface(mg_lvl%etype(sec))*(mg_lvl%esec(2,sec)-mg_lvl%esec(1,sec)+1)
         else ! 2D cells
           cvs(3) = cvs(3) + mg_lvl%esec(2,sec)-mg_lvl%esec(1,sec)+1
         endif
         call cg_npe_f(mg_lvl%etype(sec),length,ier)
         ne2vxmax=max(ne2vxmax,length)
         vx2e_size=vx2e_size+element_nvx(mg_lvl%etype(sec))*(mg_lvl%esec(2,sec)-mg_lvl%esec(1,sec)+1)
       end do
       cvs(2) = (cvs(2) + cvs(3))/2
       ne=cvs(1)
       nf=cvs(2)
       nbf=cvs(3)
       nelem=ne+nbf
       nbndry = nbf
       allocate(mg_lvl%e2vx(ne2vxmax*nelem))
       mg_lvl%nsec=nsec
       mg_lvl%ne2vx_max=ne2vxmax
       mg_lvl%nvx=nvx
       mg_lvl%nelem=nelem
       mg_lvl%nbndry=nbndry
       mg_lvl%nfaces=nf

       npmax=5
! Mesh data for cells

       call cgns_element_info(cg,mg_lvl%e2vx,mg_lvl%esec,mg_lvl%etype,mg_lvl%ne2vx,nelem,ne2vxmax,nsec)

       length = nvx
       allocate(geom%x(length))
       allocate(geom%y(length))
       allocate(geom%z(length))

       call cgns_cell_data(cg,geom%x,geom%y,geom%z,nvx)

       if (unitscale /= 1.0_8) then
         geom%x=geom%x*unitscale
         geom%y=geom%y*unitscale
         geom%z=geom%z*unitscale
       endif
!
       call mg_lvl%add_meshds(nsec,geom%x,geom%y,geom%z,vx2e_size)

! Face xyz ip coordinate and aip
       aiplength = nf*3
       allocate(geom%rip(aiplength))
       allocate(geom%aip(aiplength))

       call calc_aip_xyzip_uns(geom%x,geom%y,geom%z,geom%rip,geom%aip,mg_lvl%e2vx,mg_lvl%etype,nelem,nbndry,nvx,&
                                  ne2vxmax,mg_lvl%fine_lvl%s2g,nsec,mg_lvl%esec,nf)
! CV Center Mesh and mesh volume
       length = nelem
       allocate(geom%xc(length))
       allocate(geom%yc(length))
       allocate(geom%zc(length))

       vlength = nelem-nbndry
       allocate(geom%vol(vlength))

       call calc_vol_cv_centers_uns(geom%xc,geom%yc,geom%zc,geom%vol,geom%x,geom%y,geom%z,mg_lvl%e2vx,&
                mg_lvl%fine_lvl%gs2nb,mg_lvl%fine_lvl%gs2nb_idx,mg_lvl%etype,nvx,nelem,nbndry,nf,ne2vxmax,nsec,&
                mg_lvl%esec,geom%rip,geom%aip)

       allocate(nl(npmax))
       nl=0
       do p=1,npmax
         nl(p) = int((nelem-nbndry)/(8**p))
         if(nl(p)==0) exit
       enddo

       npmax=p-1
       print *,'max mglvl:',p

      if(agg_method(1:18)=='Oct-Tree Isotropic') call mg_lvl%generate_seeds(nl,npmax,ne,geom%vol,geom%xc,geom%yc,geom%zc)

      do p=1,npmax
        if(agg_method(1:18)=='Oct-Tree Isotropic') then
          mem_size=mg_lvl%add_transformation_bt()
        elseif(agg_method(1:20)=='Weighted Directional') then
          mem_size=mg_lvl%add_transformation(geom%xc,geom%yc,geom%zc,geom%aip,geom%vol)
        endif
        call mg_lvl%add_meshds(nsec,geom%x,geom%y,geom%z)
      enddo
    end associate

    call cg_close_f(cg%unit, ier)
    if (ier .eq. ERROR) call cg_error_exit_f

    geom%ef2nb=>geom%mg%fine_lvl%gs2nb
    geom%ef2nb_idx=>geom%mg%fine_lvl%gs2nb_idx
    geom%ne=geom%mg%fine_lvl%ng
    geom%nf=geom%mg%fine_lvl%ns
    geom%nbf=geom%mg%fine_lvl%nbs
    geom%nvx=geom%mg%nvx

   end subroutine cell_input
 end module

