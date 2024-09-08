! basic module for h3open-utilmp python api
! this module is used from h3open_py subroutines
module h3ou_py_base
  use h3ou_api, only : jcup_varp_type, jcup_varg_type
  implicit none
  private

  public :: h3ou_var_type
  public :: add_var
  public :: get_var_ptr
  public :: get_char_str
  
  integer, parameter :: STR_SHORT = 32
  integer, parameter :: STR_LONG  = 128
  integer, parameter :: STR_VLONG = 1024

  
  type h3ou_var_type
     type(jcup_varp_type), pointer :: varp_ptr
     type(jcup_varg_type), pointer :: varg_ptr
     character(len=STR_SHORT) :: var_name
     type(h3ou_var_type), pointer :: next_ptr
  end type h3ou_var_type

  type(h3ou_var_type), pointer :: varptr => null()

contains

  !=======+=========+=========+=========+=========+=========+=========+=========+

  subroutine add_var(var_name)
    implicit none
    character(len=*), intent(IN) :: var_name
    type(h3ou_var_type), pointer :: current_ptr
    
    if(.not.associated(varptr)) then
       allocate(varptr)
       varptr%varp_ptr => null()
       varptr%varg_ptr => null()
       varptr%var_name = trim(var_name)
       varptr%next_ptr => null()
       return
    end if

    current_ptr => varptr

    do while(associated(current_ptr%next_ptr))
       current_ptr => current_ptr%next_ptr
    end do

    allocate(current_ptr%next_ptr)
    current_ptr => current_ptr%next_ptr
    current_ptr%varp_ptr => null()
    current_ptr%varg_ptr => null()
    current_ptr%var_name = trim(var_name)
    current_ptr%next_ptr => null()
    
  end subroutine add_var
  
  !=======+=========+=========+=========+=========+=========+=========+=========+

  function get_var_ptr(var_name) result(res)
    implicit none
    character(len=*), intent(IN) :: var_name
    type(h3ou_var_type), pointer :: res

    res => varptr

    do while(associated(res))
       if (trim(res%var_name) == trim(var_name)) return
       res => res%next_ptr
    end do

  end function get_var_ptr
  
  !=======+=========+=========+=========+=========+=========+=========+=========+

  function get_char_str(char_str, str_len) result(res)
    implicit none
    character(len=1), intent(IN) ::char_str(*)
    integer, intent(IN) :: str_len
    character(len=STR_LONG) :: res
    integer :: i

    res = repeat(" ", STR_LONG)
    
    do i = 1, str_len
       res(i:i) = char_str(i)
    end do

  end function get_char_str
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
   
end module h3ou_py_base

  
