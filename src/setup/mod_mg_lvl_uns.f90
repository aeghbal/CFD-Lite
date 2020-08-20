
module mod_mg_lvl
    use mod_meshds
    use mod_agglomeration
    implicit none

    real, parameter :: quality_threshold=0.85

    logical , parameter , private :: lwrite_out=.false.
    integer , parameter , private :: nEigenIt = 3
    public :: mg_lvl_t
    type mg_lvl_t
        integer :: cell_no=1,nintf_c2b
        character(len=32) :: cellname
        character(len=32), allocatable, dimension(:) :: sectionName
        !private
        integer, allocatable :: esec(:,:),etype(:),sec2intf(:)
        integer, allocatable :: intf2sec(:)
        integer, allocatable :: e2vx(:),ne2vx(:)
        integer, allocatable :: vx2e(:)
        integer, allocatable :: vx2e_idx(:)
        !sec2intf is a map from section number to global interface number
        integer :: lvl_max=0,nelem,nbndry,nfaces,nvx,ne2vx_max
        type(AoP) :: gf2g(10),sf2s(10),g2gf(10),s2sf(10)
        integer :: ng_tmp=0,ns_tmp=0,nbs_tmp=0,nsec
        integer, allocatable :: gs2nb_idx_tmp(:),bs_idx_tmp(:),s2g_tmp(:)
        type(meshds_t), pointer :: fine_lvl => null()
        type(meshds_t), pointer :: last_lvl => null()
        type(meshds_t), pointer :: cur_lvl => null()
    contains
        procedure :: generate_seeds
        procedure :: map_sec2intf
        procedure :: find_element_nb
        procedure :: find_group_nb
        procedure :: goto_fine
        procedure :: goto_next
        procedure :: goto_previous
        procedure :: is_fine
        procedure :: is_last
        procedure :: add_meshds
        procedure :: add_transformation
        procedure :: add_transformation_bt
        procedure :: get_lvl
    end type mg_lvl_t

      type :: geometry_t
        integer :: ne,nf,nbf,nvx
    ! vertex coordinates
        real, pointer, dimension(:) :: x,y,z
    ! cv center coordinates
        real, pointer, dimension(:) :: xc,yc,zc
    ! integration point, area, volume
        real, pointer, dimension(:) :: rip,aip,vol
    ! element connectivity data
        type(mg_lvl_t) :: mg
        integer, pointer :: ef2nb_idx(:)=>null()
        integer, pointer :: ef2nb(:,:)=>null()
        real :: telap=0.
!        real,managed, allocatable, dimension(:) :: phic0_d,r_d,ap,b,anb,phic
!        real,allocatable, managed :: res_d(:), res_max_d(:)
!        integer, managed, allocatable :: ef2nb_d(:),ef2nb_idx_d(:)
      end type

    contains
        subroutine destroy_geom(geom)
          type(geometry_t) :: geom
          type(meshds_t), pointer :: mds
          integer :: l

          deallocate(geom%x,geom%y,geom%z,geom%xc,geom%yc,geom%zc,geom%rip,geom%aip,geom%vol)

          deallocate(geom%mg%esec)
          deallocate(geom%mg%etype)
          deallocate(geom%mg%sec2intf)
          deallocate(geom%mg%e2vx)
          deallocate(geom%mg%vx2e)
          deallocate(geom%mg%vx2e_idx)
          deallocate(geom%mg%sectionName)
          deallocate(geom%mg%intf2sec)

          call geom%mg%goto_fine()
          do
            mds => geom%mg%cur_lvl
            call destroy_meshds(mds)
            if(geom%mg%goto_next()) then
              deallocate(mds)
              l=geom%mg%cur_lvl%lvl
              deallocate(geom%mg%g2gf(l)%p,geom%mg%g2gf(l)%idx)
              deallocate(geom%mg%gf2g(l)%p)
              deallocate(geom%mg%s2sf(l)%p,geom%mg%s2sf(l)%idx)
              deallocate(geom%mg%sf2s(l)%p)
            else
              deallocate(mds)
              exit
            end if
          end do
        end subroutine

        subroutine generate_seeds(this,nseeds,nlvl,ne,vol,xc,yc,zc,phi)
          class(mg_lvl_t) :: this
          type(seed_gen_t), pointer :: seed_gen
          integer :: nlvl,nseeds(nlvl),p,ne,target_ncv,b,e
          real, dimension(ne) :: vol,xc,yc,zc
          real, pointer , optional :: phi(:)

            if(present(phi)) then
              seed_gen=>seed_gen_t(ne,vol,xc,yc,zc,phi)
            else
              seed_gen=>seed_gen_t(ne,vol,xc,yc,zc)
            endif
            do p=nlvl,1,-1
              if(p>1) then
                target_ncv = nseeds(p-1)
              else
                target_ncv = ne
              end if
              allocate(this%gf2g(p)%p(target_ncv))
              call seed_gen%grow(target_ncv,this%gf2g(p)%p,p)
            enddo
            call seed_gen%destroy_root()

        end subroutine

        subroutine add_meshds(this,nsec,x,y,z,vx2e_size)

            use mod_util
            implicit none

            class(mg_lvl_t) :: this
            integer :: s,n,ns2sf,xlen,nsec
            integer,optional:: vx2e_size
            integer :: nelem_i,nbndry_i,etype(nsec),esec(2,nsec),ne2vx_max,nvx
            integer :: nface
            type(meshds_t), pointer :: new_lvl
            character (len=32) :: cval1,cval2
            logical :: lfind
            integer :: e,idx,f,tmp,e1
            integer, allocatable :: tvec(:)
            real, dimension(:), pointer :: x,y,z

            if (present(vx2e_size)) then
                this%cell_no=1
                allocate(this%sec2intf(nsec))
                this%sec2intf=0
                allocate(this%vx2e_idx(this%nvx+1))
                allocate(this%vx2e(vx2e_size))
                nelem_i=this%nelem
                nbndry_i=this%nbndry
                nface=this%nfaces
                nsec=this%nsec
                esec=this%esec
                etype=this%etype
!
                this%lvl_max=0
                this%fine_lvl => meshds_t(0,nelem_i-nbndry_i,nface,nbndry_i,nsec)!,this%fine_lvl%previous_lvl(),this%fine_lvl%next_lvl())
                this%last_lvl => this%fine_lvl
                this%cur_lvl => this%fine_lvl
! order section ranges i.e. esec is not ordered but bs_idx is guaranteed to be ordered from low to high
                n=1
                this%fine_lvl%bs_idx(1)=1
                do n=1,nsec
                  do s=1 , nsec
                    if(this%fine_lvl%bs_idx(n)==esec(1,s)) then
                      this%fine_lvl%bs_idx(n+1)=esec(2,s)+1
                      this%etype(n) = etype(s)
                      this%esec(:,n)= esec(:,s)
                      exit
                    endif
                  enddo
                enddo
                etype = this%etype
                esec = this%esec
! boundary intefaces are actually the 2D sections of this 1 block mesh
                n=0
                do s=1,nsec
                  if(etype(s)<10) then! this is a boundary
                    n=n+1
                    this%sec2intf(s)=n
                  end if
                end do
                this%nintf_c2b=n
                allocate(this%intf2sec(n))
                do s=1,nsec
                  n=this%sec2intf(s)
                  if(this%sec2intf(s)>0) this%intf2sec(n)=s
                end do

                call this%find_element_nb(x,y,z,this%e2vx)
            else
                this%lvl_max=this%lvl_max+1
                ns2sf=this%s2sf(this%lvl_max)%idx(this%ns_tmp+1)-1
                allocate(tvec(max(this%last_lvl%ns,this%ns_tmp+1,this%last_lvl%ng,this%ng_tmp+1,ns2sf)))

                tvec(1:this%last_lvl%ns)=this%sf2s(this%lvl_max)%p(1:this%last_lvl%ns)
                deallocate(this%sf2s(this%lvl_max)%p);nullify(this%sf2s(this%lvl_max)%p)
                allocate(this%sf2s(this%lvl_max)%p(this%last_lvl%ns))
                this%sf2s(this%lvl_max)%p=tvec(1:this%last_lvl%ns)

                tvec(1:this%ns_tmp+1)=this%s2sf(this%lvl_max)%idx(1:this%ns_tmp+1)
                deallocate(this%s2sf(this%lvl_max)%idx);nullify(this%s2sf(this%lvl_max)%idx)
                allocate(this%s2sf(this%lvl_max)%idx(this%ns_tmp+1))
                this%s2sf(this%lvl_max)%idx=tvec(1:this%ns_tmp+1)
                allocate(this%s2sf(this%lvl_max)%p(ns2sf))
                this%s2sf(this%lvl_max)%p=0

                tvec(1:this%last_lvl%ng)=this%g2gf(this%lvl_max)%p(1:this%last_lvl%ng)
                deallocate(this%g2gf(this%lvl_max)%p);nullify(this%g2gf(this%lvl_max)%p)
                allocate(this%g2gf(this%lvl_max)%p(this%last_lvl%ng))
                this%g2gf(this%lvl_max)%p=tvec(1:this%last_lvl%ng)

                tvec(1:this%ng_tmp+1)=this%g2gf(this%lvl_max)%idx(1:this%ng_tmp+1)
                deallocate(this%g2gf(this%lvl_max)%idx);nullify(this%g2gf(this%lvl_max)%idx)
                allocate(this%g2gf(this%lvl_max)%idx(this%ng_tmp+1))
                this%g2gf(this%lvl_max)%idx=tvec(1:this%ng_tmp+1)

                tvec(1:this%last_lvl%ng)=this%gf2g(this%lvl_max)%p(1:this%last_lvl%ng)
                deallocate(this%gf2g(this%lvl_max)%p);nullify(this%gf2g(this%lvl_max)%p)
                allocate(this%gf2g(this%lvl_max)%p(this%last_lvl%ng))
                this%gf2g(this%lvl_max)%p=tvec(1:this%last_lvl%ng)

                new_lvl => meshds_t(this%lvl_max,this%ng_tmp,this%ns_tmp,this%nbs_tmp,this%nsec,this%last_lvl)!,new_lvl%next_lvl())

                new_lvl%s2g=this%s2g_tmp(1:this%ns_tmp)
                deallocate(this%s2g_tmp)!;nullify(this%s2g_tmp)
                new_lvl%gs2nb_idx=this%gs2nb_idx_tmp(1:this%ng_tmp+1)
                deallocate(this%gs2nb_idx_tmp)!;nullify(this%gs2nb_idx_tmp)
                new_lvl%bs_idx=this%bs_idx_tmp(1:this%nsec+1)
                deallocate(this%bs_idx_tmp)!;nullify(this%bs_idx_tmp)

                call this%last_lvl%set_next_lvl(new_lvl)
                this%last_lvl => new_lvl
                !Use ldble to detemine whethere pxyzip is single or real
                call this%find_group_nb()
            end if
            call this%last_lvl%calc_redblack()

        end subroutine add_meshds
!
        function goto_next(this) result(lsuccess)
            class(mg_lvl_t) :: this
            logical :: lsuccess

            lsuccess = associated(this%cur_lvl%next_lvl()) .and. this%cur_lvl%lvl<this%lvl_max
            if(lsuccess) this%cur_lvl => this%cur_lvl%next_lvl()

        end function goto_next

        subroutine goto_previous(this)
            class(mg_lvl_t) :: this

            if( associated(this%cur_lvl%previous_lvl())) this%cur_lvl => this%cur_lvl%previous_lvl()

        end subroutine goto_previous

        subroutine goto_fine(this)
            class(mg_lvl_t) :: this

            do while(associated(this%cur_lvl%previous_lvl()))
                this%cur_lvl => this%cur_lvl%previous_lvl()
            end do

        end subroutine goto_fine

        function is_fine(this)
            class(mg_lvl_t) :: this
            logical :: is_fine

            is_fine=.false.
            if( .not. associated(this%cur_lvl%previous_lvl())) is_fine=.true.
        end function is_fine

        function is_last(this)
            class(mg_lvl_t) :: this
            logical :: is_last

            is_last=.false.
            if( .not. associated(this%cur_lvl%next_lvl()) ) is_last=.true.
        end function is_last

        function get_lvl(this)
            class(mg_lvl_t) :: this
            integer :: get_lvl

            get_lvl=this%cur_lvl%lvl
        end function get_lvl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine find_element_nb(this,x,y,z,e2vx)
     use mod_util
     implicit none
     class(mg_lvl_t) :: this
     integer, pointer :: indices(:),f2s(:)
     integer :: s,e,j,v,e1,e2,enb,t,f,f1,f2,emax,fmax,tmp,tmp2,l,m,n,nf
     integer :: nl1,nl2,idx1,idx2,t1,t2,idx,n1,n2
     integer :: lst1(4),lst2(4),e2vx(this%ne2vx_max,this%nelem)
     real :: skewness,flatness,aspectratio
     real,dimension(*) :: x,y,z

     real , pointer :: ip(:,:)
     logical :: lfind
     logical, save :: lreported=.false.
     character (len=32) :: cval

    this%vx2e_idx=0
    this%fine_lvl%gs2nb_idx=0
    this%fine_lvl%gs2nb_idx(1)=1

    fmax=0
    do s=1,this%nsec
        do e=this%fine_lvl%bs_idx(s),this%fine_lvl%bs_idx(s+1)-1
            t=this%etype(s)
            if(e<this%fine_lvl%ng+1) this%fine_lvl%gs2nb_idx(e+1)=this%fine_lvl%gs2nb_idx(e)+element_nface(t)
            fmax=max(fmax,element_nface(t))
            do j=1,element_nvx(t)
                v=e2vx(j,e)
                this%vx2e_idx(v) = this%vx2e_idx(v) + 1
            enddo
        enddo
    enddo

    tmp=this%vx2e_idx(1)
    this%vx2e_idx(1)=1
    emax=tmp
    do v=2,this%nvx+1
        tmp2=this%vx2e_idx(v)
        emax=max(emax,tmp2)
        this%vx2e_idx(v)=this%vx2e_idx(v-1)+tmp
        tmp=tmp2
    enddo

    this%vx2e=0

    tmp=emax*fmax!maximum number of incident faces on a vertex
    allocate(indices(1:tmp))
    allocate(f2s(1:tmp))

    do s=1,this%nsec
        n1=this%fine_lvl%bs_idx(s)
        n2=this%fine_lvl%bs_idx(s+1)-1
        do e=n1,n2
            do j=1,element_nvx(this%etype(s))
                v=e2vx(j,e)
                idx1=this%vx2e_idx(v)
                do while(this%vx2e(idx1)/=0 .and. idx1<this%vx2e_idx(v+1))
                    idx1=idx1+1
                enddo
                this%vx2e(idx1)=e
            enddo
        enddo
    enddo

     nf=0!number of faces
     this%fine_lvl%gs2nb=0
     this%fine_lvl%bs=0

     do v=1,this%nvx
       indices=0
       f2s=0
       n=0
       do idx1=this%vx2e_idx(v),this%vx2e_idx(v+1)-1
         e = this%vx2e(idx1)
         do s=1,this%nsec
            if(e<this%fine_lvl%bs_idx(s+1)) exit
         enddo
         t=this%etype(s)
         do f=1,element_nface(t)
           lst1 = 0
           call face_vxlist(e2vx(:,e),t,f,lst1,nl1)
           do l=1,nl1
             if (lst1(l)==v) then
                n=n+1
                indices(n)=index_t(e,f,0)
                f2s(n)=s
             endif
           enddo
         enddo
       enddo

       do l=1,n
         do m=l+1,n
            t1=this%etype(f2s(l))
            t2=this%etype(f2s(m))
            call get_idx(indices(l),0,e1,f1)
            call get_idx(indices(m),0,e2,f2)
            if(face_compare(e1,f1,t1,e2vx(:,e1),e2,f2,t2,e2vx(:,e2))) then
                if(e1>this%fine_lvl%ng .and. e2<=this%fine_lvl%ng) then
                    idx2=this%fine_lvl%gs2nb_idx(e2)+f2-1
                    if(this%fine_lvl%gs2nb(idx2,1)/=0) cycle
                    nf=nf+1
                    this%fine_lvl%gs2nb(idx2,1)=index_t(e1,0,0)
                    this%fine_lvl%gs2nb(idx2,2)=nf
                    this%fine_lvl%s2g(nf)=indices(m)
                    this%fine_lvl%bs(e1)=index_t(e2,f2,0)
                elseif(e2>this%fine_lvl%ng .and. e1<=this%fine_lvl%ng) then
                    idx1=this%fine_lvl%gs2nb_idx(e1)+f1-1
                    if(this%fine_lvl%gs2nb(idx1,1)/=0) cycle
                    nf=nf+1
                    this%fine_lvl%gs2nb(idx1,1)=index_t(e2,0,0)
                    this%fine_lvl%gs2nb(idx1,2)=nf
                    this%fine_lvl%s2g(nf)=indices(l)
                    this%fine_lvl%bs(e2)=index_t(e1,f1,0)
                elseif(e1<=this%fine_lvl%ng .and. e2<=this%fine_lvl%ng) then
                    idx1=this%fine_lvl%gs2nb_idx(e1)+f1-1
                    idx2=this%fine_lvl%gs2nb_idx(e2)+f2-1
                    if(this%fine_lvl%gs2nb(idx1,1)/=0 .and. this%fine_lvl%gs2nb(idx2,1)/=0) cycle
                    nf=nf+1
                    this%fine_lvl%gs2nb(idx1,1)=indices(m)
                    this%fine_lvl%gs2nb(idx2,1)=indices(l)
                    this%fine_lvl%gs2nb(idx1,2)=nf
                    this%fine_lvl%s2g(nf)=indices(l)
                    if(e1>e2) then
                        this%fine_lvl%gs2nb(idx1,2)=-nf
                        this%fine_lvl%s2g(nf)=indices(m)
                    endif
                    this%fine_lvl%gs2nb(idx2,2)=-this%fine_lvl%gs2nb(idx1,2)
                else
                    if(.not. lreported) then
                      write(*,*) ' Warning: 2D element cannot be matched. '
                      lreported=.true.
                    endif
                endif
            endif
         enddo
       enddo
     enddo
!calculate edge boundary faces to prevent race condition when calcualting different interfaces
     do e1=this%fine_lvl%ng+1,this%nelem
        if(this%fine_lvl%bs(e1)<0) cycle
        call get_idx(this%fine_lvl%bs(e1),0,e,f)
        do idx=this%fine_lvl%gs2nb_idx(e),this%fine_lvl%gs2nb_idx(e+1)-1
            if(idx-this%fine_lvl%gs2nb_idx(e)+1==f) cycle
            call get_idx(this%fine_lvl%gs2nb(idx,1),0,e2,f2)
            if(f2==0) then
                this%fine_lvl%bs(e1)=-abs(this%fine_lvl%bs(e1))!because there may be at most 3 boudary faces
                this%fine_lvl%bs(e2)=-this%fine_lvl%bs(e2)
            endif
        enddo
     enddo

    allocate(ip(3,10))
    do s=1,this%nsec
        do e=this%fine_lvl%bs_idx(s),this%fine_lvl%bs_idx(s+1)-1
            t = this%etype(s)
            if(t<10) cycle
            v = element_nvx(t)

             do l = 1,v
                m = e2vx(l,e)
                ip(:,l) = [real(x(m)),real(y(m)),real(z(m))]
            enddo
            this%fine_lvl%eval(e) = esm(this%fine_lvl%evec(3*(e-1)+1:3*e),ip,v,t,skewness,flatness,aspectratio)
            !if(this%fine_lvl%eval(e)>quality_threshold) then
              !write(*,'(A,I8,A,I5,A)')         'Warning: element no. ',e,' in cell no.',this%cell_no,' has bad quality!'
              !write(*,'(A,F5.3,A,F5.3,A,F5.3)')'       skewness:',skewness,' flatness:',flatness,' aspect_ratio:',aspectratio
              !write(*,'(A,I8,I4,A,F5.3,A,F5.3,A,F5.3)') 'e:',e,t,' skewness:',skewness,' flatness:',flatness,' aspect_ratio:',aspectratio
            !end if
        enddo
    enddo
    deallocate(ip)

    if (nf/=this%fine_lvl%ns) then
      write(*,*) 'Error in creation of element neighbour list ...',nf,this%fine_lvl%ns
      do e=1,this%fine_lvl%ng
        do idx=this%fine_lvl%gs2nb_idx(e),this%fine_lvl%gs2nb_idx(e+1)-1
          if(this%fine_lvl%gs2nb(idx,2)==0) then
            do s=1,this%nsec
              if(e<this%fine_lvl%bs_idx(s+1)) exit
            enddo
            t=this%etype(s)
            f=idx-this%fine_lvl%gs2nb_idx(e)+1
            lst1 = 0
            call face_vxlist(e2vx(:,e),t,f,lst1,nl1)
            write(*,*) 'e_3D :',e,' type:',t,' sec:',s,' face:',f
            write(*,*) ' >>vxlist:',lst1(1:nl1)
          endif
        end do
      end do
      do e=this%fine_lvl%ng+1,this%fine_lvl%ng+this%fine_lvl%nbs
        if(this%fine_lvl%bs(e)==0) then
          do s=1,this%nsec
            if(e<this%fine_lvl%bs_idx(s+1)) exit
          enddo
          t=this%etype(s)
          write(*,*) 'e_2D :',e,' type:',t,' sec:',s
          write(*,*) ' >>vxlist:',e2vx(1:element_nvx(t),e)
        end if
      end do

      stop
    endif

    deallocate(f2s);deallocate(indices)
  end subroutine find_element_nb

  subroutine find_group_nb(this)
  use mod_util
  implicit none
    class(mg_lvl_t) :: this
    integer :: sec,idx,idx1,idx2,g,gf,gfnb,sf,s1,s2,g1,g2,s,n,nt
    class(meshds_t),pointer:: prv
    character(len=32) :: cval,filename
    logical :: lfind

    prv=>this%last_lvl%previous
    this%s2sf(this%lvl_max)%p=0
    do sf=1, prv%ns
        s=this%sf2s(this%lvl_max)%p(sf)
        if(s==0) cycle!because not all sf make s
        do idx=this%s2sf(this%lvl_max)%idx(abs(s)),this%s2sf(this%lvl_max)%idx(abs(s)+1)-1
            if(this%s2sf(this%lvl_max)%p(idx)==0) then
                this%s2sf(this%lvl_max)%p(idx)=sgn(s)*sf
                exit
            endif
        enddo
    enddo

    this%last_lvl%gs2nb=0
    do s=1,this%last_lvl%ns
        call get_idx(this%last_lvl%s2g(s),this%lvl_max,g1,s1)
        idx1=this%last_lvl%gs2nb_idx(g1)+s1-1
        sf=this%s2sf(this%lvl_max)%q(s,1)
        call get_idx(prv%s2g(abs(sf)),prv%lvl,gf,s2)

        if(sgn(sf)==1) then!means gf is in the g side while we seek the gf on the other side
            call get_idx(prv%gs2nb(prv%gs2nb_idx(gf)+s2-1,1),prv%lvl,gf,s2)!gf changes to be on the other side which can be 0 for boundary
        endif

        if(s2==0) then!boundary
            sec=this%nsec+1
            do while(gf<prv%bs_idx(sec))
                sec=sec-1
            enddo
            this%last_lvl%bs_idx(sec+1)=this%last_lvl%bs_idx(sec+1)-1
            this%last_lvl%gs2nb(idx1,1)=index_t(this%last_lvl%bs_idx(sec+1),0,this%lvl_max)!for bndry, gf is the number of section
            this%last_lvl%gs2nb(idx1,2)=s
            this%last_lvl%bs(this%last_lvl%bs_idx(sec+1))=index_t(g1,s1,this%lvl_max)
        else
            g2=this%gf2g(this%lvl_max)%p(gf)
            do idx2=this%last_lvl%gs2nb_idx(g2),this%last_lvl%gs2nb_idx(g2+1)-1
                if(this%last_lvl%gs2nb(idx2,1)==0) then
                    s2=idx2-this%last_lvl%gs2nb_idx(g2)+1
                    exit
                endif
            enddo
            this%last_lvl%gs2nb(idx1,1)=index_t(g2,s2,this%lvl_max)
            this%last_lvl%gs2nb(idx2,1)=index_t(g1,s1,this%lvl_max)
            this%last_lvl%gs2nb(idx1,2)=s
            if(g1>g2) this%last_lvl%gs2nb(idx1,2)=-s
            this%last_lvl%gs2nb(idx2,2)=-this%last_lvl%gs2nb(idx1,2)
        endif
    enddo
    this%last_lvl%bs_idx(1:this%nsec)=this%last_lvl%bs_idx(2:this%nsec+1)
    this%last_lvl%bs_idx(this%nsec+1)=this%last_lvl%ng+this%last_lvl%nbs+1

  end subroutine find_group_nb
!
  function add_transformation(this,xc,yc,zc,aip,vol) result(required_mem)
    use mod_util
    implicit none

    class(mg_lvl_t) :: this

    integer :: xlen,alen,vlen !Array lengths

    real, pointer :: aip(:), vol(:), xc(:), yc(:), zc(:)

    real :: gsAip(2),area,rc(3),evec(3),evec0(3),totVol,vol0,totWt,wt0,center(3)

    real , pointer :: cm(:,:),wt(:)
    integer :: sec,tmp,tmp2,s,g,gf,sf,idx,idx1,idx2,gfnb,gfnbnb,lvl,g2,l,dl,required_mem(2),ng2nb
    integer , pointer :: felem(:),ng2gf(:),gnb(:)
    integer , parameter :: g2gfmin=8,g2gfmax=12,beta=0.3
    logical :: lnew

    lvl=this%lvl_max+1

    tmp=2*this%last_lvl%ns-this%last_lvl%nbs!total number of elements in gs2nb
    allocate(wt(tmp))
    allocate(cm(3,this%last_lvl%ng))
    wt=0.0
    cm=0.
!calculate weight of different sides of a group
    do gf=1,this%last_lvl%ng
      area = 0.
      idx1 = this%last_lvl%gs2nb_idx(gf)
      idx2 = this%last_lvl%gs2nb_idx(gf+1)-1
      do idx=idx1,idx2
        s=abs(this%last_lvl%gs2nb(idx,2))
        gsAip = calc_gs_aip(aip,this%fine_lvl%ns,this%lvl_max,s,this%s2sf,this%last_lvl%evec(3*(gf-1)+1:3*gf))

        wt(idx) = gsAip(1)
        area = area + gsAip(2)
      end do
      wt(idx1:idx2)=wt(idx1:idx2)/area
      felem => fine_groups(this%g2gf,lvl-1,gf)
        do l=1,size(felem)
          cm(:,gf) = cm(:,gf) + [xc(gf),yc(gf),zc(gf)]
        enddo

      cm(:,gf)=cm(:,gf)/size(felem)
    end do

    allocate(this%gf2g(lvl)%p(1:this%last_lvl%ng))
    allocate(this%g2gf(lvl)%idx(1:this%last_lvl%ng+1))
     g=0; this%gf2g(lvl)%p = 0;this%g2gf(lvl)%idx=0

     do gf=1,this%last_lvl%ng
       if (this%gf2g(lvl)%p(gf)==0) then  ! make sure element/group is free to start a family
         g = g + 1
         totVol = 0.
         center = cm(:,gf)
         evec = this%last_lvl%evec(3*(gf-1)+1:3*gf)
         felem => fine_groups(this%g2gf,lvl-1,gf)
           do l=1,size(felem)
             totVol = totVol + vol(felem(l))
           enddo


! The parent
         this%g2gf(lvl)%idx(g)=this%g2gf(lvl)%idx(g)+1
         this%gf2g(lvl)%p(gf)  = g
!
         do idx1=this%last_lvl%gs2nb_idx(gf),this%last_lvl%gs2nb_idx(gf+1)-1
           call get_idx(this%last_lvl%gs2nb(idx1,1),this%last_lvl%lvl,gfnb,tmp)
           if (tmp==0) cycle
           rc = cm(:,gfnb) - center
           evec0= this%last_lvl%evec(3*(gfnb-1)+1:3*gfnb)
           idx2 = this%last_lvl%gs2nb_idx(gfnb)+tmp -1
           !if (this%gf2g(lvl)%p(gfnb)==0 .and. (1.+beta-abs(dot_product(rc,evec0)))>beta) then
           if (this%gf2g(lvl)%p(gfnb)==0 .and. wt(idx1)>wt(idx2) ) then

             this%g2gf(lvl)%idx(g)=this%g2gf(lvl)%idx(g)+1
             this%gf2g(lvl)%p(gfnb) = g
        ! Calcualte an estimate of the eigen vector with the addition of the new element
             felem => fine_groups(this%g2gf,lvl-1,gfnb)

             vol0=0.
               do l=1,size(felem)
                 vol0 = vol0 + vol(felem(l))
               enddo


             evec = (evec*totVol+this%last_lvl%evec(3*(gfnb-1)+1:3*gfnb)*vol0)/(totVol+vol0)
             center = (center*totVol+cm(:,gfnb)*vol0)/(totVol+vol0)
             totVol = totVol + vol0
           endif
           if (this%g2gf(lvl)%idx(g)==g2gfmin) exit
         enddo
! Add the grand children but only if the child is part of the current family
         if (this%g2gf(lvl)%idx(g)<g2gfmin) then
           do idx=this%last_lvl%gs2nb_idx(gf),this%last_lvl%gs2nb_idx(gf+1)-1
             call get_idx(this%last_lvl%gs2nb(idx,1),this%last_lvl%lvl,gfnb,tmp)
             if ( tmp/=0 .and. this%gf2g(lvl)%p(gfnb)==g) then
               do idx1=this%last_lvl%gs2nb_idx(gfnb),this%last_lvl%gs2nb_idx(gfnb+1)-1
                 call get_idx(this%last_lvl%gs2nb(idx1,1),this%last_lvl%lvl,gfnbnb,tmp)!grand kid
                 if (tmp/=0 .and. this%gf2g(lvl)%p(gfnbnb)==0) then
                   idx2 = this%last_lvl%gs2nb_idx(gfnbnb) + tmp-1
                   evec0= this%last_lvl%evec(3*(gfnbnb-1)+1:3*gfnbnb)
                   rc = (cm(:,gfnbnb) - center)
                   rc = rc/sqrt(sum(rc**2))
                   !if(abs(dot_product(evec,evec0))*(1.+beta-abs(dot_product(rc,evec0)))>beta) then
                   if (wt(idx1)>wt(idx2) ) then
                     this%g2gf(lvl)%idx(g)=this%g2gf(lvl)%idx(g)+1
                     this%gf2g(lvl)%p(gfnbnb) = g
                ! Calcualte an estimate of the eigen vector with the addition of the new element
                     felem => fine_groups(this%g2gf,lvl-1,gfnbnb)
                     vol0=0.

                       do l=1,size(felem)
                         vol0 = vol0 + vol(felem(l))
                       enddo

                     evec = (evec*totVol+evec0*vol0)/(totVol+vol0)
                     center = (center*totVol+cm(:,gfnbnb)*vol0)/(totVol+vol0)
                     totVol = totVol + vol0
                   endif
                 endif
                 if (this%g2gf(lvl)%idx(g)==g2gfmin) exit
               enddo
             endif
             if (this%g2gf(lvl)%idx(g)==g2gfmin) exit
           enddo
         endif
! At this point if there are no children attached to the parent
! add the parent to a neighboring family.
         if (this%g2gf(lvl)%idx(g)==1) then
           do idx=this%last_lvl%gs2nb_idx(gf),this%last_lvl%gs2nb_idx(gf+1)-1
             call get_idx(this%last_lvl%gs2nb(idx,1),this%last_lvl%lvl,gfnb,tmp)
             if(tmp==0) cycle
             g2 = this%gf2g(lvl)%p(gfnb)
             if(g2==0) cycle
             if (this%g2gf(lvl)%idx(g2)<g2gfmax) then
               this%gf2g(lvl)%p(gf) = g2
               this%g2gf(lvl)%idx(g2) = this%g2gf(lvl)%idx(g2)+1
               this%g2gf(lvl)%idx(g)=0
               g = g - 1
               exit
             endif
           enddo
         endif
!
       endif
     enddo

    this%ng_tmp=g
    allocate(this%g2gf(lvl)%p(1:this%last_lvl%ng));this%g2gf(lvl)%p=0

    tmp=this%g2gf(lvl)%idx(1)
    this%g2gf(lvl)%idx(1)=1
    do g=2,this%ng_tmp+1
        tmp2=this%g2gf(lvl)%idx(g)
        this%g2gf(lvl)%idx(g)=this%g2gf(lvl)%idx(g-1)+tmp
        tmp=tmp2
    enddo

    do gf=1,this%last_lvl%ng
        g=this%gf2g(lvl)%p(gf)
        do idx=this%g2gf(lvl)%idx(g),this%g2gf(lvl)%idx(g+1)-1
            if(this%g2gf(lvl)%p(idx)==0) then
                this%g2gf(lvl)%p(idx)=gf
                exit
            endif
        enddo
    enddo
    ! Calculate neighbor families (i.e. gs2nb_idx,s2sf,sf2s)
     allocate(this%sf2s(lvl)%p(1:this%last_lvl%ns))
     allocate(this%s2sf(lvl)%idx(1:this%last_lvl%ns+1))
     allocate(this%s2g_tmp(1:this%last_lvl%ns))
     allocate(this%gs2nb_idx_tmp(1:this%ng_tmp+1))
     allocate(this%bs_idx_tmp(1:this%nsec+1))
     this%bs_idx_tmp=0
     this%s2sf(lvl)%idx=0!1st just use it to count number of sides for each group

     this%sf2s(lvl)%p=0
     this%gs2nb_idx_tmp=1
     this%ns_tmp=0;this%nbs_tmp=0
!     this%ns2sf=0
     allocate(gnb(2**(min(2*lvl+num_face_bits,16))))
     do g=1,this%ng_tmp
       gnb=0
       ng2nb=this%gs2nb_idx_tmp(g)
       do idx=this%g2gf(lvl)%idx(g),this%g2gf(lvl)%idx(g+1)-1
         gf=this%g2gf(lvl)%p(idx)
         do idx2=this%last_lvl%gs2nb_idx(gf),this%last_lvl%gs2nb_idx(gf+1)-1
            call get_idx(this%last_lvl%gs2nb(idx2,1),this%last_lvl%lvl,gfnb,tmp)
            sf=this%last_lvl%gs2nb(idx2,2)
            sec=0
            if(tmp/=0) then
                gnb(ng2nb)=this%gf2g(lvl)%p(gfnb)!non boundary
            else
                do sec=1,this%nsec
                  if(gfnb<this%last_lvl%bs_idx(sec+1)) exit
                enddo
                gnb(ng2nb)=sec+this%ng_tmp!treating different boundary sections as different groups
            endif
            if ( g< gnb(ng2nb) ) then

                dl=0
                lnew=.true.
                do l=this%gs2nb_idx_tmp(g),ng2nb-1
                    if (gnb(l) == gnb(ng2nb)) then
                        lnew=.false.
                        dl=ng2nb-l-1!delta l
                        exit
                    endif
                enddo
                if(lnew) then
                    this%ns_tmp=this%ns_tmp+1
                    this%s2sf(lvl)%idx(this%ns_tmp)=this%s2sf(lvl)%idx(this%ns_tmp)+1
                    this%sf2s(lvl)%p(abs(sf))=sgn(sf)*this%ns_tmp
                    this%s2g_tmp(this%ns_tmp)=index_t(g,ng2nb,lvl)
                    if(tmp==0) then
                        this%nbs_tmp=this%nbs_tmp+1
                        this%bs_idx_tmp(sec)=this%bs_idx_tmp(sec)+1
                    else ! number of faces of adjacent group should also be increased if it is not a boundary one
                        this%gs2nb_idx_tmp(gnb(ng2nb))=this%gs2nb_idx_tmp(gnb(ng2nb))+1
                    endif
                    ng2nb=ng2nb+1
                else
                    this%sf2s(lvl)%p(abs(sf))=sgn(sf)*(this%ns_tmp-dl)
                    this%s2sf(lvl)%idx(this%ns_tmp-dl)=this%s2sf(lvl)%idx(this%ns_tmp-dl)+1
                endif
            endif
         enddo
       enddo
       this%gs2nb_idx_tmp(g)=ng2nb-1
     enddo

     deallocate(gnb)
     tmp=this%gs2nb_idx_tmp(1)
     this%gs2nb_idx_tmp(1)=1
     do g=2,this%ng_tmp+1
        tmp2=this%gs2nb_idx_tmp(g)
        this%gs2nb_idx_tmp(g)=this%gs2nb_idx_tmp(g-1)+tmp
        tmp=tmp2
     enddo
     tmp=this%s2sf(lvl)%idx(1)
     this%s2sf(lvl)%idx(1)=1
     do s=2,this%ns_tmp+1
        tmp2=this%s2sf(lvl)%idx(s)
        this%s2sf(lvl)%idx(s)=this%s2sf(lvl)%idx(s-1)+tmp
        tmp=tmp2
     enddo

     tmp=this%ng_tmp+1
     do s=2,this%nsec
        tmp2=this%bs_idx_tmp(s)
        if(this%bs_idx_tmp(s)/=0) then
            this%bs_idx_tmp(s)=this%bs_idx_tmp(s-1)+tmp
            tmp=tmp2
        endif
    enddo
    this%bs_idx_tmp(this%nsec+1)=this%bs_idx_tmp(this%nsec)+tmp2

    deallocate(wt)
! memory requirement
! snb s2g gs2nb_idx gs2nb      s2sf        s2sf_idx sf2s        g2gf        g2gf_idx gf2g
! nbs ns  ng        2(2ns-nbs) ns2sf       ns       last_lvl%ns last_lvl%ng ng       last_lvl%ng
!ns2sf=this%s2sf(lvl)%idx(this%ns_tmp+1)-1
    required_mem(1)=6*this%ns_tmp+3*this%ng_tmp+this%last_lvl%ns+this%s2sf(lvl)%idx(this%ns_tmp+1)-1+&
                2*this%last_lvl%ng-this%nbs_tmp+this%nsec+5
    required_mem(2)=4*this%ng_tmp

  end function add_transformation

  function add_transformation_bt2(this,ldble,rz,pxc,pyc,pzc,paip,pvol,xlen,alen,vlen) result(required_mem)
    use mod_util
    implicit none

    class(mg_lvl_t) :: this

    logical :: ldble !Set true if aip,vol,xc,yc,zc are in real
    real, target :: rz(*)

    integer(kind=8) :: paip, pvol, pxc, pyc, pzc
    integer :: xlen,alen,vlen !Array lengths

    real, pointer :: aip(:), vol(:), xc(:), yc(:), zc(:)

    integer :: sec,tmp,tmp2,s,g,gf,sf,idx,idx1,idx2,gfnb,gfnbnb,lvl,g2,l,dl,required_mem(2),ng2nb
    integer , pointer :: felem(:),ng2gf(:),gnb(:),gf2g_tmp(:)
    logical :: lnew

    !set up pointers based on cell precision
      aip => rz(paip:paip+alen)
      vol => rz(pvol:pvol+vlen)
      xc => rz(pxc:pxc+xlen)
      yc => rz(pyc:pyc+xlen)
      zc => rz(pzc:pzc+xlen)

    lvl=this%lvl_max+1

    allocate(this%g2gf(lvl)%p(this%last_lvl%ng))
    allocate(gf2g_tmp(this%last_lvl%ng))
    do gf=1,this%last_lvl%ng
      this%g2gf(lvl)%p(gf)=gf
      gf2g_tmp(gf)=this%gf2g(lvl)%p(gf)
    enddo

    call qsort_key(gf2g_tmp,this%g2gf(lvl)%p,1,this%last_lvl%ng)
    this%ng_tmp = gf2g_tmp(this%last_lvl%ng)! get the last coarse group no.
    allocate(this%g2gf(lvl)%idx(1:this%ng_tmp+1))
    g=0
    do idx=1,this%last_lvl%ng
      if(g/=gf2g_tmp(idx)) then
         g=gf2g_tmp(idx)
         this%g2gf(lvl)%idx(g)=idx
      end if
    end do
    this%g2gf(lvl)%idx(this%ng_tmp+1)=this%last_lvl%ng+1
    deallocate(gf2g_tmp)
    deallocate(this%gf2g(lvl)%p)
    deallocate(this%g2gf(lvl)%p)
    deallocate(this%g2gf(lvl)%idx)

  end function

  function add_transformation_bt(this) result(required_mem)
    use mod_util
    implicit none

    class(mg_lvl_t) :: this
    integer :: sec,tmp,tmp2,s,g,gf,sf,idx,idx1,idx2,gfnb,gfnbnb,lvl,g2,l,dl,required_mem(2),ng2nb
    integer , pointer :: felem(:),ng2gf(:),gnb(:),gf2g_tmp(:)
    logical :: lnew

    lvl=this%lvl_max+1
    allocate(this%g2gf(lvl)%p(this%last_lvl%ng))
    do gf=1,this%last_lvl%ng
      this%g2gf(lvl)%p(gf)=gf
    enddo
    !allocate(gf2g_tmp,source=this%gf2g(lvl)%p)
    allocate(gf2g_tmp(size(this%gf2g(lvl)%p)))
    gf2g_tmp=this%gf2g(lvl)%p

    call qsort_key_nRec(gf2g_tmp,this%g2gf(lvl)%p,this%last_lvl%ng)
    this%ng_tmp = gf2g_tmp(this%last_lvl%ng)! get the last coarse group no.
    allocate(this%g2gf(lvl)%idx(1:this%ng_tmp+1))

    g=0
    do idx=1,this%last_lvl%ng
      if(g/=gf2g_tmp(idx)) then
         g=gf2g_tmp(idx)
         this%g2gf(lvl)%idx(g)=idx
      end if
    end do
    this%g2gf(lvl)%idx(this%ng_tmp+1)=this%last_lvl%ng+1
    deallocate(gf2g_tmp)
    ! Calculate neighbor families (i.e. gs2nb_idx,s2sf,sf2s)
     allocate(this%sf2s(lvl)%p(1:this%last_lvl%ns))
     allocate(this%s2sf(lvl)%idx(1:this%last_lvl%ns+1))
     allocate(this%s2g_tmp(1:this%last_lvl%ns))
     allocate(this%gs2nb_idx_tmp(1:this%ng_tmp+1))
     allocate(this%bs_idx_tmp(1:this%nsec+1))
     this%bs_idx_tmp=0
     this%s2sf(lvl)%idx=0!1st just use it to count number of sides for each group

     this%sf2s(lvl)%p=0
     this%gs2nb_idx_tmp=1
     this%ns_tmp=0;this%nbs_tmp=0
!     this%ns2sf=0
     allocate(gnb(2**(min(2*lvl+num_face_bits,16))))
     do g=1,this%ng_tmp
       gnb=0
       ng2nb=this%gs2nb_idx_tmp(g)
       do idx=this%g2gf(lvl)%idx(g),this%g2gf(lvl)%idx(g+1)-1
         gf=this%g2gf(lvl)%p(idx)
         do idx2=this%last_lvl%gs2nb_idx(gf),this%last_lvl%gs2nb_idx(gf+1)-1
            call get_idx(this%last_lvl%gs2nb(idx2,1),this%last_lvl%lvl,gfnb,tmp)
            sf=this%last_lvl%gs2nb(idx2,2)
            sec=0
            if(tmp/=0) then
                gnb(ng2nb)=this%gf2g(lvl)%p(gfnb)!non boundary
            else
                do sec=1,this%nsec
                  if(gfnb<this%last_lvl%bs_idx(sec+1)) exit
                enddo
                gnb(ng2nb)=sec+this%ng_tmp!treating different boundary sections as different groups
            endif
            if ( g< gnb(ng2nb) ) then

                dl=0
                lnew=.true.
                do l=this%gs2nb_idx_tmp(g),ng2nb-1
                    if (gnb(l) == gnb(ng2nb)) then
                        lnew=.false.
                        dl=ng2nb-l-1!delta l
                        exit
                    endif
                enddo
                if(lnew) then
                    this%ns_tmp=this%ns_tmp+1
                    this%s2sf(lvl)%idx(this%ns_tmp)=this%s2sf(lvl)%idx(this%ns_tmp)+1
                    this%sf2s(lvl)%p(abs(sf))=sgn(sf)*this%ns_tmp
                    this%s2g_tmp(this%ns_tmp)=index_t(g,ng2nb,lvl)
                    if(tmp==0) then
                        this%nbs_tmp=this%nbs_tmp+1
                        this%bs_idx_tmp(sec)=this%bs_idx_tmp(sec)+1
                    else ! number of faces of adjacent group should also be increased if it is not a boundary one
                        this%gs2nb_idx_tmp(gnb(ng2nb))=this%gs2nb_idx_tmp(gnb(ng2nb))+1
                    endif
                    ng2nb=ng2nb+1
                else
                    this%sf2s(lvl)%p(abs(sf))=sgn(sf)*(this%ns_tmp-dl)
                    this%s2sf(lvl)%idx(this%ns_tmp-dl)=this%s2sf(lvl)%idx(this%ns_tmp-dl)+1
                endif
            endif
         enddo
       enddo
       this%gs2nb_idx_tmp(g)=ng2nb-1
     enddo

     deallocate(gnb)
     tmp=this%gs2nb_idx_tmp(1)
     this%gs2nb_idx_tmp(1)=1
     do g=2,this%ng_tmp+1
        tmp2=this%gs2nb_idx_tmp(g)
        this%gs2nb_idx_tmp(g)=this%gs2nb_idx_tmp(g-1)+tmp
        tmp=tmp2
     enddo
     tmp=this%s2sf(lvl)%idx(1)
     this%s2sf(lvl)%idx(1)=1
     do s=2,this%ns_tmp+1
        tmp2=this%s2sf(lvl)%idx(s)
        this%s2sf(lvl)%idx(s)=this%s2sf(lvl)%idx(s-1)+tmp
        tmp=tmp2
     enddo

     tmp=this%ng_tmp+1
     do s=2,this%nsec
        tmp2=this%bs_idx_tmp(s)
        if(this%bs_idx_tmp(s)/=0) then
            this%bs_idx_tmp(s)=this%bs_idx_tmp(s-1)+tmp
            tmp=tmp2
        endif
    enddo
    this%bs_idx_tmp(this%nsec+1)=this%bs_idx_tmp(this%nsec)+tmp2
! memory requirement
! snb s2g gs2nb_idx gs2nb      s2sf        s2sf_idx sf2s        g2gf        g2gf_idx gf2g
! nbs ns  ng        2(2ns-nbs) ns2sf       ns       last_lvl%ns last_lvl%ng ng       last_lvl%ng
!ns2sf=this%s2sf(lvl)%idx(this%ns_tmp+1)-1
    required_mem(1)=6*this%ns_tmp+3*this%ng_tmp+this%last_lvl%ns+this%s2sf(lvl)%idx(this%ns_tmp+1)-1+&
                2*this%last_lvl%ng-this%nbs_tmp+this%nsec+5
    required_mem(2)=4*this%ng_tmp

  end function

  subroutine map_sec2intf(this,sec,intf)
    class(mg_lvl_t) :: this
    integer :: sec,intf
    this%sec2intf(sec)=intf
  end subroutine map_sec2intf

  recursive function calc_gs_aip(aip,ns,lvl,s,s2sf,vec) result(gs_aip)
    real :: gs_aip(2)
    integer , intent(in) :: ns,lvl,s
    real , intent(in) :: aip(3,ns)
    real :: vec(3),area,norm(3)
    integer :: idx
    type(AoP) , intent(in) :: s2sf(10)

    gs_aip=0.

    if(lvl>0) then
        do idx=s2sf(lvl)%idx(s),s2sf(lvl)%idx(s+1)-1
            gs_aip=gs_aip+calc_gs_aip(aip,ns,lvl-1,abs(s2sf(lvl)%p(idx)),s2sf,vec)!*sgn(s2sf(idx,lvl))
        end do
    else
        gs_aip(2)=sqrt(aip(1,s)**2+aip(2,s)**2+aip(3,s)**2)
        norm = aip(:,s)/gs_aip(2)
        gs_aip(1)=1.0-abs(norm(1)*vec(1)+norm(2)*vec(2)+norm(3)*vec(3))
    endif
  end function calc_gs_aip

  recursive subroutine n_fine_sides(s2sf,lvl,s,n)
    integer :: n
    integer , intent(in) :: lvl,s
    type (AoP) , intent(in) :: s2sf(10)
    integer :: idx

    if(lvl > 0) then
        do idx=s2sf(lvl)%idx(s),s2sf(lvl)%idx(s+1)-1
            call n_fine_sides(s2sf,lvl-1,abs(s2sf(lvl)%p(idx)),n)
        end do
    else
        n=n+1
    endif
  end subroutine n_fine_sides

  recursive subroutine retrieve_fsides(s2sf,lvl,s,fsides,sz,cnt,side_sgn)
  use mod_util
  implicit none
    integer , intent(in) :: lvl,s,sz
    integer ,intent(inout):: fsides(sz),cnt,side_sgn
    type (AoP) , intent(in) :: s2sf(10)
    integer :: idx

    if(lvl > 1) then
        do idx=s2sf(lvl)%idx(s),s2sf(lvl)%idx(s+1)-1
            side_sgn=sgn(s2sf(lvl)%p(idx))*side_sgn
            call retrieve_fsides(s2sf,lvl-1,abs(s2sf(lvl)%p(idx)),fsides,sz,cnt,side_sgn)
        end do
    else
        do idx=s2sf(1)%idx(s),s2sf(1)%idx(s+1)-1
            fsides(cnt)=side_sgn*s2sf(1)%p(idx)
            cnt=cnt+1
        enddo
    endif
  end subroutine retrieve_fsides


  function fine_sides(s2sf,lvl,s) result(fsides)
  use mod_util
  implicit none

  integer :: sz , lvl , s , cnt,s_sgn
  type(AoP) :: s2sf(10)
  integer , pointer :: fsides(:)

  sz=0
  call n_fine_sides(s2sf,lvl,s,sz)
  allocate(fsides(1:sz))

  cnt=1
  s_sgn=sgn(s)
  call retrieve_fsides(s2sf,lvl,s,fsides,sz,cnt,s_sgn)

  end function fine_sides

  recursive subroutine n_fine_groups(g2gf,lvl,g,n)
    integer , intent(inout):: n
    integer , intent(in) :: lvl,g
    type (AoP) , intent(in) :: g2gf(10)
    integer :: idx

    if(lvl > 0) then
        do idx=g2gf(lvl)%idx(g),g2gf(lvl)%idx(g+1)-1
            call n_fine_groups(g2gf,lvl-1,g2gf(lvl)%p(idx),n)
        end do
    else
        n=n+1
    endif
  end subroutine n_fine_groups

  recursive subroutine retrieve_fgroups(g2gf,lvl,g,fgroups,sz,cnt)
    integer , intent(in) :: lvl,g,sz
    integer ,intent(inout):: fgroups(sz),cnt
    type (AoP) , intent(in) :: g2gf(10)
    integer :: idx
!here lvl is transformation level from 1 to 10
    if(lvl > 0) then
        do idx=g2gf(lvl)%idx(g),g2gf(lvl)%idx(g+1)-1
            call retrieve_fgroups(g2gf,lvl-1,g2gf(lvl)%p(idx),fgroups,sz,cnt)
        end do
    else
        fgroups(cnt) = g
        cnt = cnt +1
    endif
  end subroutine retrieve_fgroups

  function fine_groups(g2gf,lvl,g) result(fgroups)

  integer :: sz , lvl , g , cnt
  type(AoP) :: g2gf(10)
  integer , pointer :: fgroups(:)

  sz=0
  call n_fine_groups(g2gf,lvl,g,sz)
  allocate(fgroups(1:sz))

  cnt=1
  call retrieve_fgroups(g2gf,lvl,g,fgroups,sz,cnt)

  end function fine_groups

  function coarse_side(sf2s,lvl,sf) result(s)
  use mod_util
  implicit none
    integer :: sf, lvl , l , s
    type(AoP) :: sf2s(10)

    s=sf
    do l=1,lvl
        s=sgn(s)*sf2s(l)%p(abs(s))
    enddo
 end function coarse_side

 function coarse_group(gf2g,lvl,gf) result(g)
    integer :: gf, lvl , l , g
    type(AoP) :: gf2g(10)

    g=gf
    do l=1,lvl
        g=gf2g(l)%p(abs(g))
    enddo
 end function coarse_group
! Equiangle skewness measure(esm)
 function esm(evec,vx,nvx,t,skewness,flatness,aspectRatio) result(eval)
    use mod_util
    implicit none
    real :: skewness,aspectRatio,flatness, eval,evec(3)
    real ,parameter :: pi=3.14
    integer :: nvx,t,v
    real :: vx(3,nvx),angle(4),a(3,4),a_min(3),a_max(3),cv_center(3),cs_center(3),CF,CVX,face_norm(3),ab(3),ac(3)
    real :: d(4),d_min,d_max,theta_min,theta_max,theta0!equiangle
    integer :: vxlist(10),lst(4),f,n,i,imin,imax

    vxlist = [1,2,3,4,5,6,7,8,9,10]
    flatness = 0.
    skewness = 0.
    aspectRatio =0.
    d_min = 1000.
    d_max = 0.
    cv_center = sum(vx,2)/element_nvx(t)

    if(t==TETRA_4 .or. t==PYRA_5 .or. t==PENTA_6 .or. t==HEXA_8) then
     do f=1,element_nface(t)
      lst = 0
      call face_vxlist(vxlist,t,f,lst,n)
      if(n==3) then
        a(:,1)=vx(:,lst(2))-vx(:,lst(1));d(1)=sqrt(sum(a(:,1)**2));a(:,1) = a(:,1)/d(1)
        a(:,2)=vx(:,lst(3))-vx(:,lst(2));d(2)=sqrt(sum(a(:,2)**2));a(:,2) = a(:,2)/d(2)
        a(:,3)=vx(:,lst(1))-vx(:,lst(3));d(3)=sqrt(sum(a(:,3)**2));a(:,3) = a(:,3)/d(3)
        angle(1)=acos(dot_product(-a(:,1),a(:,2)))*180./pi
        angle(2)=acos(dot_product(-a(:,2),a(:,3)))*180./pi
        angle(3)=acos(dot_product(-a(:,3),a(:,1)))*180./pi
        theta0 = 60
        cs_center = (vx(:,lst(1))+vx(:,lst(2))+vx(:,lst(3)))/3.
      elseif(n==4) then
        a(:,1)=vx(:,lst(2))-vx(:,lst(1));d(1)=sqrt(sum(a(:,1)**2));a(:,1) = a(:,1)/d(1)
        a(:,2)=vx(:,lst(3))-vx(:,lst(2));d(2)=sqrt(sum(a(:,2)**2));a(:,2) = a(:,2)/d(2)
        a(:,3)=vx(:,lst(4))-vx(:,lst(3));d(3)=sqrt(sum(a(:,3)**2));a(:,3) = a(:,3)/d(3)
        a(:,4)=vx(:,lst(1))-vx(:,lst(4));d(4)=sqrt(sum(a(:,4)**2));a(:,4) = a(:,4)/d(4)
        angle(1)=acos(dot_product(-a(:,1),a(:,2)))*180./pi
        angle(2)=acos(dot_product(-a(:,2),a(:,3)))*180./pi
        angle(3)=acos(dot_product(-a(:,3),a(:,4)))*180./pi
        angle(4)=acos(dot_product(-a(:,4),a(:,1)))*180./pi
        theta0 = 90
        cs_center = (vx(:,lst(1))+vx(:,lst(2))+vx(:,lst(3))+vx(:,lst(4)))/4.
      endif
      theta_max = maxval(angle(1:n))
      theta_min = minval(angle(1:n))
      skewness = max(skewness,max((theta_max-theta0)/(180-theta0),(theta0-theta_min)/theta0))
! calculate flatness of the face(centroid test)
      CVX=0.
      do v=1,n
        CVX = max(CVX,sqrt(dot_product(cv_center-vx(:,lst(v)),cv_center-vx(:,lst(v)))))
      enddo
      ab=vx(:,lst(2))-vx(:,lst(1))
      ac=vx(:,lst(3))-vx(:,lst(1))
      call vec_a_bcrossc(face_norm,ab,ac,1)! outward
      face_norm = face_norm/sqrt(dot_product(face_norm,face_norm))
      ab =cs_center-cv_center
      ab = ab/sqrt(dot_product(ab,ab))
      flatness = max(flatness,1.-dot_product(ab,face_norm))
! calculate Aspect Ratio
      imin=n+1
      do i=1,n
        if(d(i)<d_min) imin=i
      enddo
      if(imin<=n) then
        d_min=d(imin)
        a_min=a(:,imin)
      endif

      imax=n+1
      do i=1,n
        if(d(i)>d_max) imax=i
      enddo
      if(imax<=n) then
        d_max=d(imax)
        a_max=a(:,imax)
      endif
     enddo
! Find the vertical component of d_min to d_max to calculate aspect ratio
     a(:,1)=a_min-a_max*dot_product(a_min,a_max)
     d_min=d_min*dot_product(a(:,1),a_min)
     aspectRatio = 1.-d_min/d_max
    endif

    eval = max(skewness,aspectRatio,flatness)

    evec = [maxval(vx(1,:))-minval(vx(1,:)),maxval(vx(2,:))-minval(vx(2,:)),maxval(vx(3,:))-minval(vx(3,:))]
    evec = evec/sqrt(sum(evec**2))
  end function
!
! Principal Component Analysis
 function pca(evec,vx,nvx) result(eval)
    real :: eval , evec(3)
    integer :: nvx
    real :: vx(3,nvx)
    real :: cov(3,3)

    integer :: v,i,j,k
!    logical :: ldegen(3)
    real :: s(3),tr

    s = sum(vx,2)/nvx
    do v = 1, nvx
      vx(:,v) = vx(:,v)-s
    enddo

! create covariance tensor
    cov = 0.
    do i=1,3
      do j=1,3
        if(i<=j) then
          do v=1,nvx
            cov(i,j)=cov(i,j)+vx(i,v)*vx(j,v)
          enddo
        else
          cov(i,j) = cov(j,i)
        endif
      enddo
    enddo
    cov = cov / (nvx-1)

    tr = 0.! trace of covariance tensor
    do i =1,3
      tr = tr + cov(i,i)
    enddo
! Power Iteration method to find principal eigenvector
    evec = [1.,1.,1.]/sqrt(3.)
    do k = 1, nEigenIt
      s = 0.
      do i = 1,3
        do j = 1,3
          s(i) = s(i)+cov(i,j)*evec(j)
        enddo
      enddo
      evec = s/sqrt(sum(s**2))
    enddo
! calculate eigenvalue

!    ldegen = .false.
! degeneracy check
! x-y isotropy ldegen(1)=T
! x-z isotropy ldegen(2)=T
! y-z isotropy ldegen(3)=T
    s = 0
    do i = 1,3
      do j = 1,3
	     s(i) = s(i)+cov(i,j)*evec(j)
!	     if(i<j .and. abs(evec(i) - evec(j))<1e-5) ldegen(i+j-2) = .true.
      enddo
    enddo

    eval = 0.
	 do k = 1,3
      eval = eval + evec(k)*s(k)
    enddo

    evec = abs(evec)
! Normalize eigen value
    eval = eval / tr
  end function


end module mod_mg_lvl

