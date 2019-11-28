module mod_eqn_setup
  use mod_cell
  use mod_properties
  implicit none


  ! Boundary conditions
  type :: bc_t
    character(len=32) :: bc_type
    integer :: idx! interface number
    integer :: esec(2)! range of halo elements for this bc
    character(len=32) :: name ! bc name
    procedure (i_intfSub) , pointer :: coef => null()
  end type

  type :: equation_t
    ! name
    character(len=16) :: name
    ! Fields
    real, pointer, dimension(:) :: phi,phi0
    ! gradient
    real, pointer, dimension(:) :: grad
    ! boundary conditions
    type(bc_t), allocatable, dimension(:) :: bcs

   contains
    procedure :: make_bc
  end type

  interface
    subroutine i_intfSub(bc,geom,eqn,prop)
      import bc_t
      import geometry_t
      import equation_t
      import properties_t

      class(bc_t) :: bc
      type(geometry_t) :: geom
      class(equation_t) :: eqn
      type(properties_t) :: prop

    end subroutine i_intfSub
  end interface

  contains
  function make_bc(eqn,geom,secName) result(bc)
    class(equation_t),target :: eqn
    type(bc_t), pointer :: bc
    character(len=*) :: secName
    type(geometry_t) :: geom
    integer :: i,s

    nullify(bc)
    do i=1,geom%mg%nintf_c2b
      s=geom%mg%intf2sec(i)
      if(index(trim(secName),trim(geom%mg%sectionName(s)))>0) then
        bc=>eqn%bcs(i)
        bc%idx=i
        bc%name=secName
        bc%esec=geom%mg%esec(:,s)
      end if
    end do
    if(.not. associated(bc)) then
      print *,'BC '//trim(secName)//' not found!'
      stop
    end if

  end function

end module
