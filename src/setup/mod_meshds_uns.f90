module mod_meshds
    implicit none

    type AoP!Array of Pointers
        integer , pointer :: p(:)=>null()
        integer , pointer :: idx(:)=>null()
      contains
        procedure :: q
    end type AoP
    interface assignment (=)
       module procedure assign_aop
    end interface

    type meshds_t
        !private
        integer :: lvl=0,ns=0,ng=0,nbs=0,nred=0!number of global sides(ngs) and number of groups(ng)
        integer , pointer :: gs2nb(:,:)=>null(),gs2nb_idx(:)=>null(),s2g(:)=>null(),bs(:)=>null(),bs_idx(:)=>null()
        integer , pointer :: red(:)=>null(),black(:)=>null()
        real , pointer :: eval(:)=>null(),evec(:)=>null()!eigenValue,eigenVector
        type(meshds_t),pointer :: next=>null(),previous=>null()
      contains
        procedure :: next_lvl
        procedure :: previous_lvl
        procedure :: set_next_lvl
        procedure :: calc_redblack
    end type meshds_t
    interface meshds_t
        module procedure constructor_meshds_t
    end interface meshds_t

    interface assignment(=)
        module procedure assign_meshds
    end interface

    contains

        function constructor_meshds_t(level,ng,ns,nbs,nsec,prv,nxt) result(meshds)
    !prv and nxt are for not putting target when calling them which is wrong
    !this should create transformations
            type(meshds_t),pointer :: meshds
            type(meshds_t),pointer,optional :: nxt,prv
            integer :: level,ng,ns,nbs,nsec,nred

            allocate(meshds)

            meshds%ng=ng
            meshds%ns=ns
            meshds%nbs=nbs
            meshds%nred=ng

            allocate(meshds%gs2nb_idx(ng+1))
            allocate(meshds%bs_idx(nsec+1))
            allocate(meshds%gs2nb(1:ns*2-nbs,1:2))
            allocate(meshds%s2g(ns+1))
            allocate(meshds%bs(ng+1:ng+nbs))
            allocate(meshds%red(ng+1))
            allocate(meshds%eval(ng))
            allocate(meshds%evec(3*ng))

			if(present(prv)) meshds%previous=>prv
			if(present(nxt)) meshds%next=>nxt
			meshds%lvl=level

        end function constructor_meshds_t

        subroutine destroy_meshds(meshds)
            type(meshds_t), pointer :: meshds

            deallocate(meshds%gs2nb_idx)
            deallocate(meshds%bs_idx)
            deallocate(meshds%gs2nb)
            deallocate(meshds%s2g)
            deallocate(meshds%bs)
            deallocate(meshds%red)
            deallocate(meshds%eval)
            deallocate(meshds%evec)

        end subroutine

        subroutine assign_meshds(mdso,mdsi)
            type(meshds_t) , intent(inout) :: mdso
            type(meshds_t) , intent(in) :: mdsi

            mdso%lvl=mdsi%lvl
            mdso%ns=mdsi%ns
            mdso%ng=mdsi%ng
            mdso%nbs=mdsi%nbs
            mdso%nred=mdsi%nred
            mdso%gs2nb=>mdsi%gs2nb
            mdso%gs2nb_idx=>mdsi%gs2nb_idx
            mdso%s2g=>mdsi%s2g
            mdso%bs=>mdsi%bs
            mdso%bs_idx=>mdsi%bs_idx
            mdso%red=>mdsi%red
            mdso%black=>mdsi%black
            mdso%eval=>mdsi%eval
            mdso%evec=>mdsi%evec
        end subroutine

        function next_lvl(this)
            class(meshds_t) :: this
            type(meshds_t),pointer :: next_lvl

            next_lvl=>this%next
        end function

        function previous_lvl(this)
            class(meshds_t) :: this
            type(meshds_t),pointer :: previous_lvl

            previous_lvl=>this%previous
        end function previous_lvl

        subroutine set_next_lvl(this,nxt)
            class(meshds_t) :: this
            type(meshds_t) , pointer :: nxt

            this%next=>nxt
        end subroutine set_next_lvl

        subroutine assign_aop(a,b)!a=(/g,s,val/)
            type(aop) , intent(inout) :: a
            integer , intent(in) :: b(3)
            a%p(a%idx(b(1))+b(2)-1)=b(3)
        end subroutine assign_aop

        function q(this,g,s) result(res)
            class(aop) :: this
            integer :: res
            integer :: g,s
            res=g/abs(g)*this%p(this%idx(abs(g))+s-1)
        end function q

        subroutine calc_redblack(this)
            use mod_util
            implicit none
            class(meshds_t) :: this
            integer :: g,gnb,idx,fl,r,b
            integer, pointer :: rb(:)
            logical :: lfound,lexit

            this%nred=0
            this%red=0
            allocate(rb(1:this%ng))
            rb=0

            g=1
            lexit=.false.
            do while(g<=this%ng)
    !finding an element that is not either red/black
                do while(rb(g)/=0)
                    if(g==this%ng) then
                        lexit=.true.
                        exit
                    endif
                    g=g+1
                enddo
                if(lexit) exit
    !figuring out if neighbors of it are not black
                lfound=.true.
                do idx=this%gs2nb_idx(g),this%gs2nb_idx(g+1)-1
                    call get_idx(this%gs2nb(idx,1),this%lvl,gnb,fl)
                    if(fl/=0) then
                        if(rb(gnb)==2) then
                            lfound=.false.
                            exit
                        endif
                    endif
                enddo
    !if all were red or boundary and not black save neighbors as red and center as black
                if(lfound) then
                    do idx=this%gs2nb_idx(g),this%gs2nb_idx(g+1)-1
                        call get_idx(this%gs2nb(idx,1),this%lvl,gnb,fl)
                        if(fl/=0 ) then
                            if(rb(gnb) == 0) then
                                this%nred=this%nred+1
                                rb(gnb)=1
                            endif
                        endif
                    enddo
                    rb(g)=2
                else
                    this%nred=this%nred+1
                    rb(g)=1
                endif
                g=g+1
            enddo

            b=this%nred
            this%nred=0
            do g=1,this%ng
                if(rb(g)==1) then!red
                    this%nred=this%nred+1
                    this%red(this%nred)=g
                elseif(rb(g)==2) then !black
                    b=b+1
                    this%red(b)=g
                endif
            enddo
            deallocate(rb)
            nullify(this%black)
            this%black=>this%red(this%nred+1:this%ng)
            this%red=>this%red(1:this%nred)
        end subroutine calc_redblack


end module mod_meshds
