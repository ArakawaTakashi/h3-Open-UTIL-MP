module mod_comp
  use jcup_interface, only : jcup_varp_type, jcup_varg_type
  use mod_utils, only : STR_SHORT, STR_MID, STR_LONG
  implicit none
  private
  
!--------------------------------   public  ----------------------------------!

  public :: varp_type
  public :: get_num_of_varp          ! integer function (self)
  public :: get_varp_ptr             ! type(varp_type), pointer function (self, varp_name or varp_num)
  public :: varg_type
  public :: get_num_of_varg          ! integer function (self)
  public :: get_varg_ptr             ! type(varg_type), pointer function (self, varg_name or varg_num)
  public :: comp_type
  public :: init_comp                ! subroutine (comp_name)
  public :: get_num_of_comp          ! integer function ()
  public :: get_comp_ptr             ! type(comp_type), pointer function (comp_name or comp_num)
  
!--------------------------------   private  ---------------------------------!

  interface get_varp_ptr
     module procedure get_varp_ptr_name, get_varp_ptr_num
  end interface get_varp_ptr

  interface get_varg_ptr
     module procedure get_varg_ptr_name, get_varg_ptr_num
  end interface get_varg_ptr

  interface get_comp_ptr
     module procedure get_comp_ptr_name, get_comp_ptr_num
  end interface get_comp_ptr
  
  type varp_type ! put data type
     type(jcup_varp_type), pointer :: jcup_varp_ptr
     character(len=STR_SHORT) :: name
     character(len=STR_SHORT) :: grid_name
     integer                  :: num_of_layer = 1
     type(varp_type), pointer :: next_ptr 
  end type varp_type

  type varg_type ! get data type
     type(jcup_varg_type), pointer :: jcup_varg_ptr
     character(len=STR_SHORT) :: name
     character(len=STR_SHORT) :: my_grid_name
     character(len=STR_SHORT) :: send_comp
     character(len=STR_SHORT) :: send_grid
     character(len=STR_SHORT) :: send_data
     integer                  :: time_lag = -1     ! 0 or -1 or 1
     integer                  :: data_tag = 1      ! data tag is used for grid interpolation 
     integer                  :: intvl = 3600
     integer                  :: num_of_layer = 1
     character(len=3)         :: flag ="SNP"       ! "SNP" or "AVR"
     type(varg_type), pointer :: next_ptr 
  end type varg_type

  type comp_type
     character(len=STR_SHORT) :: name
     integer :: num_of_put_data = 0
     type(varp_type), pointer :: pdata_ptr => null()
     integer :: num_of_get_data = 0
     type(varg_type), pointer :: gdata_ptr => null()
     type(comp_type), pointer :: next_ptr
  end type comp_type

  integer :: num_of_comp = 0
  type(comp_type), pointer :: comp_ptr => null()
  
contains

  !=======+=========+=========+=========+=========+=========+=========+=========+
  subroutine add_varp(self, name, grid_name, num_of_layer, num_of_data)
    implicit none
    type(varp_type), pointer :: self
    character(len=*), intent(IN) :: name
    character(len=*), intent(IN) :: grid_name
    integer, intent(IN)          :: num_of_layer
    integer, intent(INOUT)       :: num_of_data
    type(varp_type), pointer :: bptr
    type(varp_type), pointer :: cptr

    cptr => self
    bptr => null()


    do while (associated(cptr))
       if (trim(cptr%name) == trim(name)) return
       bptr => cptr
       cptr => cptr%next_ptr
    end do

    allocate(cptr)

    if (associated(bptr)) then
       bptr%next_ptr => cptr
    else
       self => cptr
    end if

    cptr%next_ptr => null()
    cptr%name         = trim(name)
    cptr%grid_name    = trim(grid_name)
    cptr%num_of_layer = num_of_layer
    num_of_data       = num_of_data + 1
    
  end subroutine add_varp
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_varp_ptr_name(self, name) result(cptr)
    implicit none
    type(varp_type), pointer :: self
    character(len=*), intent(IN) :: name
    type(varp_type), pointer :: cptr

    cptr => self

    do while(associated(cptr))
       if (trim(cptr%name) == trim(name)) return
       cptr => cptr%next_ptr
    end do

  end function get_varp_ptr_name
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_varp_ptr_num(self, varp_num) result(cptr)
    implicit none
    type(varp_type), pointer :: self
    integer, intent(IN)      :: varp_num
    type(varp_type), pointer :: cptr
    integer :: counter
    
    cptr => self
    counter = 1
    
    do while(associated(cptr))
       if (varp_num == counter) return
       cptr => cptr%next_ptr
       counter = counter + 1
    end do

  end function get_varp_ptr_num
  
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  subroutine add_varg(self, name, grid_name, send_comp, send_grid, send_data, num_of_layer, time_lag, &
                      data_tag, intvl, flag, num_of_data)
    implicit none
    type(varg_type), pointer :: self
    character(len=*), intent(IN) :: name
    character(len=*), intent(IN) :: grid_name
    character(len=*), intent(IN) :: send_comp
    character(len=*), intent(IN) :: send_grid
    character(len=*), intent(IN) :: send_data
    integer, intent(IN)          :: num_of_layer
    integer, intent(IN)          :: time_lag
    integer, intent(IN)          :: data_tag
    integer, intent(IN)          :: intvl
    character(len=3), intent(IN) :: flag
    integer, intent(INOUT)       :: num_of_data
    type(varg_type), pointer :: bptr
    type(varg_type), pointer :: cptr

    cptr => self
    bptr => null()

    do while (associated(cptr))
       bptr => cptr
       cptr => cptr%next_ptr
    end do

    allocate(cptr)

    if (associated(bptr)) then
       bptr%next_ptr => cptr
    else
       self => cptr
    end if

    cptr%next_ptr => null()
    cptr%name          = trim(name)
    cptr%my_grid_name  = trim(grid_name)
    cptr%send_comp     = trim(send_comp)
    cptr%send_grid     = trim(send_grid)
    cptr%send_data     = trim(send_data)
    cptr%num_of_layer  = num_of_layer
    cptr%time_lag      = time_lag
    cptr%data_tag      = data_tag
    cptr%flag          = flag
    cptr%intvl         = intvl

    num_of_data = num_of_data + 1
    
  end subroutine add_varg
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_varg_ptr_name(self, name) result(cptr)
    implicit none
    type(varg_type), pointer :: self
    character(len=*), intent(IN) :: name
    type(varg_type), pointer :: cptr

    cptr => self

    do while(associated(cptr))
       if (trim(cptr%name) == trim(name)) return
       cptr => cptr%next_ptr
    end do

  end function get_varg_ptr_name
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_varg_ptr_num(self, varg_num) result(cptr)
    implicit none
    type(varg_type), pointer :: self
    integer, intent(IN)      :: varg_num
    type(varg_type), pointer :: cptr
    integer :: counter
    
    cptr => self
    counter = 1
    
    do while(associated(cptr))
       if (varg_num == counter) return
       cptr => cptr%next_ptr
       counter = counter + 1
    end do

  end function get_varg_ptr_num
  
  
  !=======+=========+=========+=========+=========+=========+=========+=========+

  subroutine init_comp(name)
    implicit none
    character(len=*), intent(IN) :: name
    type(comp_type), pointer :: cptr

    call add_comp(comp_ptr, name)
    
  end subroutine init_comp
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  subroutine add_comp(self, name)
    implicit none
    type(comp_type), pointer :: self
    character(len=*), intent(IN) :: name
    type(comp_type), pointer :: bptr
    type(comp_type), pointer :: cptr

    cptr => self
    bptr => null()


    do while (associated(cptr))
       if (trim(cptr%name) == trim(name)) return
       bptr => cptr
       cptr => cptr%next_ptr
    end do

    allocate(cptr)

    if (associated(bptr)) then
       bptr%next_ptr => cptr
    else
       self => cptr
    end if

    num_of_comp = num_of_comp + 1
    cptr%next_ptr => null()
    cptr%name = trim(name)
    call set_comp(cptr, name)

  end subroutine add_comp
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_num_of_comp() result(res)
    implicit none
    integer :: res

    res = num_of_comp

  end function get_num_of_comp

  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_comp_ptr_name(name) result(cptr)
    implicit none
    character(len=*), intent(IN) :: name
    type(comp_type), pointer :: cptr

    cptr => comp_ptr

    do while(associated(cptr))
       if (trim(cptr%name) == trim(name)) return
       cptr => cptr%next_ptr
    end do

  end function get_comp_ptr_name
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_comp_ptr_num(comp_num) result(cptr)
    implicit none
    integer, intent(IN)      :: comp_num
    type(comp_type), pointer :: cptr
    integer :: counter
    
    cptr => comp_ptr
    counter = 1
    
    do while(associated(cptr))
       if (comp_num == counter) return
       cptr => cptr%next_ptr
       counter = counter + 1
    end do

  end function get_comp_ptr_num
  
  !=======+=========+=========+=========+=========+=========+=========+=========+

  subroutine set_comp(self, comp_name)
    use mod_utils, only : put_log, put_error
    use mod_namelist, only : GRID_LATLON, GRID_ICO, &
                             get_num_of_configuration, &
                             get_put_comp_name, get_put_grid_name, &
                             get_get_comp_name, get_get_grid_name, &
                             exchange_data_type, &
                             get_num_of_put_data, get_ed_ptr_from_put_comp_and_num, &
                             get_num_of_recv_target, &
                             get_num_of_get_data, get_ed_ptr_from_get_comp_and_num
    implicit none
    type(comp_type), pointer :: self
    character(len=*), intent(IN) :: comp_name
    type(exchange_data_type), pointer :: ed_ptr
    integer :: i, j
    character(len=STR_LONG) :: log_str
    type(varp_type), pointer :: vpptr
    type(varg_type), pointer :: vgptr

    self%name = trim(comp_name)

    self%num_of_put_data = 0

    do i = 1, get_num_of_put_data(trim(comp_name))
       ed_ptr => get_ed_ptr_from_put_comp_and_num(trim(comp_name), i)
       call add_varp(self%pdata_ptr, trim(ed_ptr%var_put), trim(ed_ptr%parent%put_grid_name), &
                     ed_ptr%num_of_layer, &
                     self%num_of_put_data)
    end do
    
    self%num_of_get_data = 0
    
    do i = 1, get_num_of_recv_target(trim(comp_name))
       do j = 1, get_num_of_get_data(trim(comp_name), i)
          ed_ptr => get_ed_ptr_from_get_comp_and_num(trim(comp_name), i, j)
          call add_varg(self%gdata_ptr, trim(ed_ptr%var_get), trim(ed_ptr%parent%get_grid_name), &
                        trim(ed_ptr%parent%put_comp_name), trim(ed_ptr%parent%put_grid_name), &
                        trim(ed_ptr%var_put), ed_ptr%num_of_layer, ed_ptr%lag, ed_ptr%grid_intpl_tag, ed_ptr%intvl, &
                        ed_ptr%flag, &
                        self%num_of_get_data)
       end do
    end do

    call put_log("")
    call put_log("--------------- component information")
    call put_log("component name : "//trim(comp_name))

    write(log_str, '(A,I6)') "num of put data = ", self%num_of_put_data
    call put_log(trim(log_str))

    do i = 1, self%num_of_put_data
       vpptr => get_varp_ptr(self%pdata_ptr, i)
       write(log_str, '(A,A10,A,A10)') "    name : ",trim(vpptr%name),", grid : ", trim(vpptr%grid_name)
       call put_log(trim(log_str))
    end do

    write(log_str, '(A,I6)') "num of get data = ", self%num_of_get_data
    call put_log(trim(log_str))

    do i = 1, self%num_of_get_data
       vgptr => get_varg_ptr(self%gdata_ptr, i)
       write(log_str, '(A,A10,A,A10)') "    name : ",trim(vgptr%name),", grid : ", trim(vgptr%my_grid_name)
       call put_log(trim(log_str))
       write(log_str, '(A,A)')  "        send comp : ",trim(vgptr%send_comp)
       call put_log(trim(log_str))
       write(log_str, '(A,A)')  "        send grid : ",trim(vgptr%send_grid)
       call put_log(trim(log_str))
       write(log_str, '(A,A)')  "        send data : ",trim(vgptr%send_data)
       call put_log(trim(log_str))
       write(log_str, '(A,I6)') "        layer     : ",vgptr%num_of_layer
       call put_log(trim(log_str))
       write(log_str, '(A,I6)') "        data tag  : ",vgptr%data_tag
       call put_log(trim(log_str))
       write(log_str, '(A,I6)') "        interval  : ",vgptr%intvl
       call put_log(trim(log_str))
    end do
       
  end subroutine set_comp
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_num_of_varp(self) result(res)
    implicit none
    type(comp_type), pointer :: self
    integer :: res

    res = self%num_of_put_data

  end function get_num_of_varp
  
  !=======+=========+=========+=========+=========+=========+=========+=========+
  function get_num_of_varg(self) result(res)
    implicit none
    type(comp_type), pointer :: self
    integer :: res

    res = self%num_of_get_data

  end function get_num_of_varg
  

end module mod_comp
