module mod_agglomeration
  use mod_dll
  implicit none

  type :: bb_t
    integer :: parent_group_id=0,group_id=-1
    real :: rmin(3)
    real :: rmax(3)
    real :: group_vol ! Sum V in this bb / target_ncv

    type(dll_t) :: list
  endtype

 type :: btree_node_t
   logical :: isLeftThread=.true.
   logical :: isRightThread=.true.
   integer :: id=0
   type(bb_t), pointer :: val=>null()
!   type(btree_node_t), pointer :: parent=>null()
   type(btree_node_t), pointer :: left=>null()
   type(btree_node_t), pointer :: right=>null()
 end type

 type :: tbt_iterator_t
   type(btree_node_t), pointer :: node=>null()
  contains
   procedure :: left_turn
   procedure :: right_turn
   procedure :: isLeaf
 end type

 interface assignment(=)
   module procedure assign_val_btree
 end interface

 interface btree_node_t
   module procedure construct_btree_node_t
 end interface

  type :: seed_gen_t
    integer :: ncv=0,nleaf=0
    integer :: target_ncv=1
    logical :: lphi = .false.
    real, pointer :: vol(:)=>null()
    real, pointer :: phi(:)=>null()
    real, pointer :: xc(:)=>null()
    real, pointer :: yc(:)=>null()
    real, pointer :: zc(:)=>null()
    real :: vol_ave
    integer :: last_id=0
    type(btree_node_t), pointer :: root=>null()
    type(btree_node_t), pointer :: begin=>null()
    type(btree_node_t), pointer :: end=>null()
  contains
    procedure :: rewind_tbt
    procedure :: split_node
    procedure :: iterate_next
    procedure :: iterate_prev
    procedure :: goto_next_leaf
    procedure :: destroy_root
    procedure :: destroy_leaf
    procedure :: goto_1st_leaf
    procedure :: destroy_root_seed
    procedure :: destroy_leaf_seed
    procedure :: grow
    procedure :: split_leaf
    procedure :: get_seeds
  end type


  interface seed_gen_t
    procedure construct_seed_gen
  end interface

 contains
 function construct_btree_node_t(id,val) result(node)
  type(btree_node_t), pointer :: node
  type(bb_t), pointer :: val
  integer :: id

  allocate(node)
  node%id = id
  node%val => val

 end function

 subroutine assign_val_btree(val,iterator)
  type(bb_t), pointer, intent(out) :: val
  type(tbt_iterator_t), intent(in) :: iterator

  val => iterator%node%val
 end subroutine

function left_turn(this) result(lsuccess)
   class(tbt_iterator_t) :: this
   logical :: lsuccess

   lsuccess = .not. this%node%isLeftThread
   if(lsuccess) this%node => this%node%left
 end function

 function right_turn(this) result(lsuccess)
   class(tbt_iterator_t) :: this
   logical :: lsuccess

   lsuccess = .not. this%node%isRightThread
   if(lsuccess) this%node => this%node%right
 end function

 logical function isLeaf(this)
   class(tbt_iterator_t) :: this

   isLeaf = this%node%isLeftThread .and. this%node%isRightThread

 end function

!!!!!!!!!!!!!!!!!!!!!!!
! tree member functions
!!!!!!!!!!!!!!!!!1!!!!!
 function rewind_tbt(this,iterator) result(ldone)
   class(seed_gen_t) :: this
   type(tbt_iterator_t) :: iterator
   logical :: ldone

   ldone = associated(this%root)
   if(ldone) iterator%node => this%begin

 end function

 subroutine split_node(this,it,left_val,right_val)
  class(seed_gen_t) :: this
  type(tbt_iterator_t) :: it
  type(bb_t),optional, pointer :: left_val,right_val
  type(btree_node_t), pointer :: left=>null(),right=>null()

  if(.not. (it%node%isLeftThread .and. it%node%isRightThread) ) then
    print *,'This node is not a leaf.'
!    pause
    return
  end if

  left => it%node%left
  right => it%node%right

  if(present(left_val)) then
    it%node%isLeftThread = .false.
    this%last_id = this%last_id +1
    it%node%left => btree_node_t(this%last_id,left_val)
    it%node%left%left => left
    it%node%left%right => it%node
    if(associated(it%node,target=this%begin)) this%begin => it%node%left
  endif
  if(present(right_val)) then
    it%node%isRightThread = .false.
    this%last_id = this%last_id +1
    it%node%right=> btree_node_t(this%last_id,right_val)
    it%node%right%left => it%node
    it%node%right%right => right
    if(associated(it%node,target=this%end)) this%end => it%node%right
  endif
  if( .not. present(left_val) .and. .not. present(right_val)) then
    print *,'Split error!'
    stop
  endif
 end subroutine

 function iterate_next(this,iterator) result(lsuccess)
   class(seed_gen_t) :: this
   type(tbt_iterator_t) :: iterator
   logical :: lsuccess

   lsuccess= .not. associated(iterator%node,target=this%end)
   if(.not. lsuccess) return
   if(iterator%right_turn()) then
     do while(iterator%left_turn())
     end do
   else
     iterator%node => iterator%node%right
   endif

 end function

  function iterate_prev(this,iterator) result(lsuccess)
   class(seed_gen_t) :: this
   type(tbt_iterator_t) :: iterator
   logical :: lsuccess

   lsuccess= .not. associated(iterator%node,target=this%begin)
   if(.not. lsuccess) return
   if(iterator%left_turn()) then
     do while(iterator%right_turn())
     end do
   else
     iterator%node => iterator%node%left
   endif

 end function

 function goto_next_leaf(this,iterator) result(lsuccess)
    class(seed_gen_t) :: this
    type(tbt_iterator_t) :: iterator
    logical :: lsuccess


    lsuccess=.false.
    do while(this%iterate_next(iterator))
      lsuccess=iterator%isLeaf()
      if(lsuccess) exit
    end do

 end function

 subroutine destroy_root(this)
    class(seed_gen_t) :: this
    type(tbt_iterator_t) :: iterator
    logical :: ldestroyed,lsuccess


    if(this%rewind_tbt(iterator)) then
      do while(associated(iterator%node))
        ldestroyed = this%destroy_leaf(iterator)
        if( .not. ldestroyed ) lsuccess= this%iterate_next(iterator)
      end do
    end if

 end subroutine

 function destroy_leaf(this,it) result(lsuccess)
  class(seed_gen_t) :: this
  type(tbt_iterator_t) :: it
  type(btree_node_t), pointer :: left=>null(),right=>null()
  logical :: lsuccess

  lsuccess = it%node%isLeftThread .and. it%node%isRightThread
  if(.not. lsuccess) return

  left => it%node%left
  right => it%node%right

  if(it%node%id==1) then
    ! root node
    deallocate(it%node)
    nullify(it%node)
  elseif(associated(it%node,target=right%left)) then
    ! left leaf removal
    if(associated(it%node,target=this%begin)) this%begin => right
    deallocate(it%node)
    it%node=>right
    it%node%left => left
    it%node%isLeftThread = .true.
  elseif(associated(it%node,target=left%right)) then
    ! right leaf removal
    if(associated(it%node,target=this%end)) this%end => left
    deallocate(it%node)
    it%node=>left
    it%node%right=>right
    it%node%isRightThread = .true.
  endif

 end function

  function goto_1st_leaf(this,it) result(lsuccess)
    class(seed_gen_t) :: this
    type(tbt_iterator_t) :: it
    logical :: lsuccess

    nullify(it%node)
    lsuccess =this%rewind_tbt(it)
    if(lsuccess) then
      if(.not. isLeaf(it)) lsuccess = this%goto_next_leaf(it)
    end if
  end function

  subroutine destroy_root_seed(this)
    class(seed_gen_t) :: this
    type(tbt_iterator_t) :: iterator
    logical :: ldestroyed,lsuccess

    if(this%goto_1st_leaf(iterator)) then
      do while(associated(iterator%node))
        ldestroyed = this%destroy_leaf_seed(iterator)
        if( .not. ldestroyed ) lsuccess= this%iterate_next(iterator)
      end do
    end if
    if(.not. this%lphi) then
      deallocate(this%phi)
    end if

  end subroutine

 function destroy_leaf_seed(this,it) result(lsuccess)
  class(seed_gen_t) :: this
  type(tbt_iterator_t) :: it
  type(btree_node_t), pointer :: left=>null(),right=>null()
  class(bb_t), pointer :: bb
  logical :: lsuccess

  lsuccess = it%node%isLeftThread .and. it%node%isRightThread
  if(.not. lsuccess) return

  left => it%node%left
  right => it%node%right

  if(it%node%id==1) then
    ! root node
    deallocate(it%node)
    nullify(it%node)
  elseif(associated(it%node,target=right%left)) then
    ! left leaf removal
    if(associated(it%node,target=this%begin)) this%begin => right
    bb = it
    call bb%list%destroy()
    deallocate(it%node)
    it%node=>right
    it%node%left => left
    it%node%isLeftThread = .true.
  elseif(associated(it%node,target=left%right)) then
    ! right leaf removal
    if(associated(it%node,target=this%end)) this%end => left
    bb = it
    call bb%list%destroy()
    deallocate(it%node)
    it%node=>left
    it%node%right=>right
    it%node%isRightThread = .true.
  endif

 end function


  function construct_seed_gen(ncv,vol,xc,yc,zc,phi) result(seed_gen)
    type(seed_gen_t), pointer :: seed_gen
    integer :: ncv,target_ncv
    real, target, dimension(ncv) :: vol,xc,yc,zc
    real, pointer, optional :: phi(:)
    type(bb_t), pointer :: root_bb
    class(value_t), pointer :: val
    real :: d
    integer :: e

    allocate(root_bb)
    allocate(seed_gen)
    if(present(phi)) then
      seed_gen%phi=>phi
      seed_gen%lphi = .true.
    else
      allocate(seed_gen%phi(ncv))
      seed_gen%phi = 1.
    end if
    seed_gen%last_id = 1
    seed_gen%nleaf=1
    seed_gen%root => btree_node_t(seed_gen%last_id,root_bb)
    seed_gen%root%left => seed_gen%root
    seed_gen%root%right=> seed_gen%root
    seed_gen%begin => seed_gen%root
    seed_gen%end => seed_gen%root

    seed_gen%ncv=ncv
    seed_gen%target_ncv=1
    seed_gen%vol=>vol
    seed_gen%xc=>xc
    seed_gen%yc=>yc
    seed_gen%zc=>zc
    ! add root bb
    root_bb%rmin = [minval(xc),minval(yc),minval(zc)]
    root_bb%rmax = [maxval(xc),maxval(yc),maxval(zc)]

    seed_gen%vol_ave = 0.
    do e=1,ncv
      val=>value_t(e)
      call root_bb%list%push(val)
      seed_gen%vol_ave = seed_gen%vol_ave + vol(e)*seed_gen%phi(e)
    end do

    root_bb%group_vol = 2.1*seed_gen%vol_ave!1.
    root_bb%group_id=1

  end function

  function split_leaf(this,tbt_it,threshold_density,nchild) result(lsuccess)
    class(seed_gen_t) :: this
    logical :: lsuccess
    integer :: nchild
    type(tbt_iterator_t) :: tbt_it
    type(iterator_t) :: it
    type(bb_t), pointer :: bb
    type(bb_t), pointer :: left,right
    class(value_t), pointer :: val
    integer :: i,isplit,e
    real, pointer :: pos(:)=>null()
    real, dimension(3) :: lrmin,rrmin,lrmax,rrmax
    real :: d,vol_l,vol_r,threshold_density

    nchild = 0
    lsuccess=.false.
    if(.not. isLeaf(tbt_it)) return
    bb = tbt_it
    ! density = bb%group_vol/this%vol_ave
    if(bb%group_vol/this%vol_ave<threshold_density) return
    if(bb%list%length==1) return
    lsuccess = .true.! split when density greater than threshold denisity. Typically when it is greater than 2
    ! split in the longest direction
    d=0.
    isplit = 0
    do i=1,3
      if(abs(bb%rmin(i)-bb%rmax(i))>d)then
        d=abs(bb%rmin(i)-bb%rmax(i))
        isplit = i
      end if
    end do
    d = (bb%rmin(isplit)+bb%rmax(isplit))/2. ! find the middle of the bb
    allocate(left)
    allocate(right)
    left%rmin=bb%rmin
    left%rmax=bb%rmax
    right%rmin=bb%rmin
    right%rmax=bb%rmax

    select case(isplit)
     case(1)
       pos => this%xc
       left%rmax(1)=d
       right%rmin(1)=d
     case(2)
       pos => this%yc
       left%rmax(2)=d
       right%rmin(2)=d
     case(3)
       pos => this%zc
       left%rmax(3)=d
       right%rmin(3)=d
    end select

    lrmin=1e20
    lrmax=-1e20
    rrmin=1e20
    rrmax=-1e20

    vol_l=0.
    vol_r=0.
    if(bb%list%rewind(it)) then
      do while(bb%list%pop(it,val))
        e=val%id
        if(pos(e)<d) then
          call left%list%push(val)
          vol_l=vol_l+this%vol(e)*this%phi(e)
          lrmin(1)=min(lrmin(1),this%xc(e))
          lrmin(2)=min(lrmin(2),this%yc(e))
          lrmin(3)=min(lrmin(3),this%zc(e))
          lrmax(1)=max(lrmax(1),this%xc(e))
          lrmax(2)=max(lrmax(2),this%yc(e))
          lrmax(3)=max(lrmax(3),this%zc(e))
        else
          call right%list%push(val)
          vol_r=vol_r+this%vol(e)*this%phi(e)
          rrmin(1)=min(rrmin(1),this%xc(e))
          rrmin(2)=min(rrmin(2),this%yc(e))
          rrmin(3)=min(rrmin(3),this%zc(e))
          rrmax(1)=max(rrmax(1),this%xc(e))
          rrmax(2)=max(rrmax(2),this%yc(e))
          rrmax(3)=max(rrmax(3),this%zc(e))
        end if
      end do
    end if

    left%group_vol = vol_l
    left%rmin=lrmin
    left%rmax=lrmax

    right%group_vol = vol_r
    right%rmin=rrmin
    right%rmax=rrmax

    if(bb%group_id==-1) then! this is an intermediate node
      left%parent_group_id=bb%parent_group_id
      right%parent_group_id=bb%parent_group_id
    else
      left%parent_group_id=bb%group_id
      right%parent_group_id=bb%group_id
    endif

    if(left%list%length>0 .and. right%list%length>0) then
      nchild = 2
      call this%split_node(tbt_it,left_val=left,right_val=right)
    elseif(left%list%length>0) then
      nchild = 1
      call this%split_node(tbt_it,left_val=left)
      deallocate(right)
    elseif(right%list%length>0) then
      nchild = 1
      call this%split_node(tbt_it,right_val=right)
      deallocate(left)
    else
      print *,"Error: This parent node does not have leaves!"
      stop
    endif

  end function

  subroutine grow(this,new_target_ncv,gf2g,lvl)
    class(seed_gen_t) :: this
    integer :: new_target_ncv
    type(tbt_iterator_t) :: it
    type(iterator_t) :: list_it
    type(bb_t), pointer :: leaf
    integer :: gf2g(new_target_ncv)
    class(value_t), pointer :: val
    real :: threshold_density
    integer :: nleaf,nleaf0,nchild,lvl
    logical :: lsuccess

    nleaf0= this%nleaf
    threshold_density = 2.
    this%vol_ave = this%vol_ave*(this%target_ncv)/(new_target_ncv)
    this%target_ncv = new_target_ncv

    if(lvl>1) then
     do while (this%nleaf<new_target_ncv)
      if(this%goto_1st_leaf(it)) then
      do
        if(this%split_leaf(it,threshold_density,nchild)) then
          lsuccess= this%goto_next_leaf(it)
          this%nleaf = this%nleaf+nchild-1
          if(this%nleaf==new_target_ncv) exit
        endif
        if(.not. this%goto_next_leaf(it)) exit
      end do
      end if
      if(this%nleaf==nleaf0) threshold_density=threshold_density*0.75
      nleaf0=this%nleaf
     enddo
    ! now number the leafs which are the groups in this lvl
     nullify(it%node)
     nleaf=0
     if(this%goto_1st_leaf(it)) then
      do
        nleaf = nleaf + 1
        leaf = it
        ! if a leaf does not split, its element number will still change
        if(leaf%group_id/=-1) leaf%parent_group_id=leaf%group_id
        leaf%group_id = nleaf
        gf2g(nleaf) = leaf%parent_group_id
        if(.not. this%goto_next_leaf(it)) exit
      end do
     endif
    elseif(lvl==1) then
    ! now number the leafs which are the groups in this lvl
     if(this%goto_1st_leaf(it)) then
      do
        leaf = it
        if(leaf%list%rewind(list_it)) then
          do while(leaf%list%pop(list_it,val))
            gf2g(val%id) = leaf%group_id
            deallocate(val)
          end do
        endif
        if(.not. this%goto_next_leaf(it)) exit
      end do
     endif
    endif
  end subroutine

  subroutine get_seeds(this,seeds)
    class(seed_gen_t) :: this
    integer :: seeds(this%target_ncv)
    type(tbt_iterator_t) :: tbt_it
    type( bb_t), pointer :: leaf
    type(iterator_t) :: it
    type(value_t),pointer :: val
    real :: d,dmin,dr(3),center(3)
    integer :: n,e,emin

    ! pick the cv that is closer to the center of the leaf
    n=0
    if(this%goto_1st_leaf(tbt_it)) then
    do
      !call assign_val_btree(leaf,tbt_it)
      leaf => tbt_it%node%val
      center = (leaf%rmin+leaf%rmax)*0.5
      dmin=1e10
      if(leaf%list%rewind(it)) then
      do
        !val = it
        val => it%node%val
        e=val%id
        dr=[this%xc(e)-center(1),this%yc(e)-center(2),this%zc(e)-center(3)]
        d=sqrt(dot_product(dr,dr))

        if(dmin>d) then
          dmin=d
          emin=e
        end if
        if(.not. leaf%list%iterate(it)) exit
      end do
      end if
      n=n+1
      seeds(n) = emin
      if(.not. this%goto_next_leaf(tbt_it)) exit
    enddo
    endif

  end subroutine

end module

