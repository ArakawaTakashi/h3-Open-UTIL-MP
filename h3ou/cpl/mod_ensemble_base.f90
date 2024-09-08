!====================================================================================================
!> @brief
!> h3open_mp ensemble base modele
!
module mod_ensemble_base
  use mpi
  implicit none
  private

!--------------------------------   public  ----------------------------------!
  integer, parameter, public :: NO_ENSEMBLE = 0
  integer, parameter, public :: ONE_TO_ONE  = 1
  integer, parameter, public :: MANY_TO_ONE = 2
  
  integer, parameter, public :: STR_LEN = 1024
  integer, parameter, public :: NAME_LEN = 64

  type mpi_info_type
     character(len=NAME_LEN) :: name
     integer :: comm          ! mpi communicator
     integer :: comm_id       ! communicator id number
     integer :: leader_rank   ! global rank of my leader
     integer :: rank
     integer :: group
     integer :: size
  end type


  type(mpi_info_type), pointer, public :: global => null()
  type(mpi_info_type), pointer, public :: ecomp  => null() 
  type(mpi_info_type), pointer, public :: glocal => null()
  type(mpi_info_type), pointer, public :: local  => null()
  type(mpi_info_type), pointer, public :: couple => null()

  public :: eb_init                 ! subroutine (comp_name, num_of_ensemble)
  public :: get_ensemble_type       ! integer function()
  public :: set_coupling_flag       ! subroutine (coupling_flag)
  public :: get_coupling_flag       ! logical function 
  public :: set_ecomp_rank          ! subroutine (comp_name, comp_id, start_rank, next_rank)
  public :: set_local_name          ! subroutine ()
  public :: set_target_name         ! subroutine ()
  public :: get_my_name             ! character(len=NAME_LEN) function ()
  public :: get_target_name         ! character(len=NAME_LEN) function ()
  public :: get_couple_id           ! integer function ()
  public :: eb_end_init             ! subroutine ()
  public :: get_local_mpi_param     ! subroutine (comm, group, size, rank)
  public :: get_local_mpi_rank      ! integer function ()
  public :: is_local_leader         ! logical function ()
  public :: get_ensemble_mean       ! subroutine (data, is_get_ok)
  public :: bcast_ensemble          ! subroutine (data)  
  public :: send_local              ! subroutine
  public :: recv_local              ! subroutine
  
!--------------------------------   private  ---------------------------------!

  interface get_ensemble_mean
     module procedure get_ensemble_mean_1d, get_ensemble_mean_2d
  end interface get_ensemble_mean

  interface bcast_ensemble
     module procedure bcast_ensemble_int
     module procedure bcast_ensemble_1d, bcast_ensemble_2d
  end interface bcast_ensemble

  interface send_local
     module procedure send_int_1d_local, send_real_1d_local, send_double_1d_local
     module procedure send_int_2d_local, send_real_2d_local, send_double_2d_local
     module procedure send_int_3d_local, send_real_3d_local, send_double_3d_local
  end interface send_local

  interface recv_local
     module procedure recv_int_1d_local, recv_real_1d_local, recv_double_1d_local
     module procedure recv_int_2d_local, recv_real_2d_local, recv_double_2d_local
     module procedure recv_int_3d_local, recv_real_3d_local, recv_double_3d_local
  end interface recv_local
  
  integer, parameter, private :: MPI_MY_TAG = 0

  integer :: ensemble_type    ! NO_ENSEMBLE, ONE_TO_ONE, ONE_TO_MANY

  logical :: coupling_flag = .true. ! flag of this component is coupled or not
  
  character(len=NAME_LEN) :: my_name          ! my component name
  character(len=NAME_LEN) :: target_name      ! coupled component name
  
  integer :: ierror
  integer :: errorcode

  real(kind=8), pointer :: recv_buffer_1d(:)
  real(kind=8), pointer :: recv_buffer_2d(:)
  
contains

!=======+=========+=========+=========+=========+=========+=========+=========+
!> @breaf
!> initialize module
subroutine eb_init(comp_name, num_of_ensemble)
  implicit none
  character(len=*), intent(IN) :: comp_name 
  integer, optional, intent(IN) :: num_of_ensemble
  character(len=NAME_LEN) :: e_comp_name ! ensemble component name
  logical :: init_flag
  integer :: start_rank, next_rank
  integer :: glocal_id
  integer :: i
  
  call MPI_initialized(init_flag, ierror)

  if (.not.init_flag) call MPI_init(ierror)
  
  allocate(global)
  
  global%name = "GLOBAL"
  global%comm_id = 1
  global%comm = MPI_COMM_WORLD
  global%leader_rank = 0
  
  call mpi_comm_size(global%comm, global%size, ierror)
  call mpi_comm_rank(global%comm, global%rank, ierror)
  call mpi_comm_group(global%comm, global%group, ierror)
  
  call set_ensemble_type(num_of_ensemble)

  allocate(ecomp)
  allocate(glocal)
  allocate(couple)
  allocate(local)

end subroutine eb_init

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_ensemble_type(num_of_ensemble)
  implicit none
  integer, intent(IN) :: num_of_ensemble
  integer :: int_array(1), recv_array(1)
  integer :: min, max

  int_array(1) = num_of_ensemble

  call mpi_allreduce(int_array, recv_array, 1, MPI_INTEGER, MPI_MIN, global%comm, ierror)
  min = recv_array(1)

  call mpi_allreduce(int_array, recv_array, 1, MPI_INTEGER, MPI_MAX, global%comm, ierror)
  max = recv_array(1)

  if ((min == 1).and.(max == 1)) then
     ensemble_type = NO_ENSEMBLE
  else if ((min == 1).and.(max > 1)) then
     ensemble_type = MANY_TO_ONE
  else if ((min > 1).and.(max > 1)) then
     ensemble_type = ONE_TO_ONE
    if (min /= max) then
       write(0, *) "ERROR!!!, ensemble number mismatch"
       call MPI_abort(global%comm, errorcode, ierror)
    end if
 end if
 
end subroutine set_ensemble_type


!=======+=========+=========+=========+=========+=========+=========+=========+

function get_ensemble_type() result(res)
  implicit none
  integer :: res

  res = ensemble_type

end function get_ensemble_type

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_coupling_flag(cflag)
  implicit none
  logical, intent(IN) :: cflag

  coupling_flag = cflag

end subroutine set_coupling_flag

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_coupling_flag() result(res)
  implicit none
  logical :: res

  res = coupling_flag

end function get_coupling_flag


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_ecomp_rank(comp_name, ecomm_id, start_rank, next_rank)
  implicit none
  character(len=*), intent(IN) :: comp_name
  integer, intent(IN) :: ecomm_id
  integer, intent(IN) :: start_rank
  integer, intent(OUT) :: next_rank
  character(len=NAME_LEN) :: name_buffer
  integer :: int_array(1), recv_array(1)

  if (global%rank == start_rank) then
     int_array(1) = len_trim(comp_name)
     name_buffer = trim(comp_name)
  end if
 
  call mpi_bcast(int_array, 1, MPI_INTEGER, start_rank, global%comm, ierror)
  call mpi_bcast(name_buffer, int_array(1), MPI_CHARACTER, start_rank, global%comm, ierror)

  name_buffer = name_buffer(1:int_array(1))
  
  if (trim(name_buffer) == trim(ecomp%name)) then
     int_array(1) = 1
  else
     int_array(1) = 0
  end if

  call mpi_allreduce(int_array, recv_array, 1, MPI_INTEGER, MPI_SUM, global%comm, ierror)

  next_rank = recv_array(1) + start_rank

  if ((global%rank >= start_rank).and.(global%rank < next_rank)) then
     ecomp%comm_id = ecomm_id
     ecomp%leader_rank = start_rank
     ecomp%size = next_rank - start_rank
     ecomp%rank = global%rank - start_rank
  end if

end subroutine set_ecomp_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_local_name(num_of_ensemble)
  implicit none
  integer, intent(IN) :: num_of_ensemble
  integer :: my_ensemble_id
  

  local%size = ecomp%size/num_of_ensemble
  
  my_ensemble_id = ecomp%rank/local%size

  couple%comm_id = my_ensemble_id + 1
  local%comm_id  = (ecomp%comm_id-1)*num_of_ensemble + couple%comm_id
  
  if (ensemble_type == NO_ENSEMBLE) then
     local%name = trim(ecomp%name)
  else
    write(local%name, '(A,I0.4)') trim(ecomp%name),my_ensemble_id
 end if
 
  my_name = trim(local%name)
  
end subroutine set_local_name

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_target_name()
  implicit none
  integer :: int_array(1), recv_array(1)
  character(len=NAME_LEN) :: name_buffer

  if (.not.get_coupling_flag()) then
     target_name = "NO_TARGET"
     return
  end if
  
  if (couple%rank == 0) then
     int_array(1) = len_trim(my_name)
     name_buffer = trim(my_name)
  end if
 
  call mpi_bcast(int_array, 1, MPI_INTEGER, 0, couple%comm, ierror)
  call mpi_bcast(name_buffer, int_array(1), MPI_CHARACTER, 0, couple%comm, ierror)

  name_buffer = name_buffer(1:int_array(1))

  if (trim(name_buffer) /= trim(my_name)) target_name = trim(name_buffer)

  if (couple%rank == couple%size-1) then
     int_array(1) = len_trim(my_name)
     name_buffer = trim(my_name)
  end if
 
  call mpi_bcast(int_array, 1, MPI_INTEGER, couple%size-1, couple%comm, ierror)
  call mpi_bcast(name_buffer, int_array(1), MPI_CHARACTER, couple%size-1, couple%comm, ierror)

  name_buffer = name_buffer(1:int_array(1))

  if (trim(name_buffer) /= trim(my_name)) target_name = trim(name_buffer)

end subroutine set_target_name

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine eb_end_init(is_debug_mode)
  use jcup_interface, only : jcup_set_world, jcup_set_new_comp, jcup_initialize, &
                             jcup_get_mpi_parameter, jcup_log
  implicit none
  logical, intent(IN) :: is_debug_mode
  integer :: log_level
  integer :: i
  
  call jcup_set_world(couple%comm)
  
  call jcup_set_new_comp(trim(local%name))

  if (is_debug_mode) then
     log_level = 2
  else
     log_level = 0
  end if

  
  call jcup_initialize(trim(local%name), "SEC", log_level, .false.)

  call jcup_get_mpi_parameter(trim(local%name), local%comm, local%group, local%size, local%rank)

  glocal%comm_id = (local%rank + 1)*10 + ecomp%comm_id

  call mpi_comm_split(global%comm, glocal%comm_id, global%rank, glocal%comm, ierror)
  call mpi_comm_size(glocal%comm, glocal%size, ierror)
  call mpi_comm_rank(glocal%comm, glocal%rank, ierror)
  call mpi_comm_group(glocal%comm, glocal%group, ierror)

  my_name = trim(local%name)
  
  !do i = 1, jcup_get_num_of_component()
  !   target_name = jcup_get_component_name(i)
  !   if (trim(target_name) /= trim(my_name)) exit
  !end do
    
end subroutine eb_end_init

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_my_name() result(res)
  implicit none
  character(len=NAME_LEN) :: res

  res = trim(my_name)

end function get_my_name

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_target_name() result(res)
  implicit none
  character(len=NAME_LEN) :: res

  res = trim(target_name)

end function get_target_name

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_couple_id() result(res)
  implicit none
  integer :: res

  res = couple%comm_id

end function get_couple_id

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_local_mpi_param(comm, group, size, rank)
  implicit none
  integer, intent(OUT) :: comm, group, size, rank

  comm  = local%comm
  group = local%group
  size  = local%size
  rank  = local%rank

end subroutine get_local_mpi_param

 
!=======+=========+=========+=========+=========+=========+=========+=========+

function get_local_mpi_rank() result(res)
  implicit none
  integer :: res

  res = local%rank

end function get_local_mpi_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

function is_local_leader() result(res)
  implicit none
  logical :: res

  res = (local%rank == 0)

end function is_local_leader

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_ensemble_mean_1d(data, is_get_ok)
  implicit none
  real(kind=8), intent(INOUT) :: data(:)
  logical, intent(OUT) :: is_get_ok
  integer :: data_size
  
  data_size = size(data)

  if (.not.associated(recv_buffer_1d)) then
     allocate(recv_buffer_1d(data_size))
  else if (size(recv_buffer_1d) < data_size) then
     deallocate(recv_buffer_1d)
     allocate(recv_buffer_1d(data_size))
  end if
  
  call MPI_Reduce(data, recv_buffer_1d, data_size, MPI_DOUBLE, MPI_SUM, 0, glocal%comm, ierror)

  if (glocal%rank == 0) then
     data(:) = recv_buffer_1d(:)/glocal%size
     is_get_ok = .true.
  else
     is_get_ok = .false.
  end if
   
end subroutine get_ensemble_mean_1d

  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_ensemble_mean_2d(data, is_get_ok)
  implicit none
  real(kind=8), intent(INOUT) :: data(:,:)
  logical, intent(OUT) :: is_get_ok
  integer :: size1, size2
  integer :: i, j, counter
  
  size1 = size(data,1)
  size2 = size(data,2)

  if (.not.associated(recv_buffer_2d)) then
     allocate(recv_buffer_2d(size1*size2))
  else if (size(recv_buffer_2d) < size1*size2) then
     deallocate(recv_buffer_2d)
     allocate(recv_buffer_2d(size1*size2))
  end if
  
  call MPI_Reduce(data, recv_buffer_2d, size1*size2, MPI_DOUBLE, MPI_SUM, 0, glocal%comm, ierror)

  if (glocal%rank == 0) then
     counter = 0
     do j = 1, size2
        do i = 1, size1
           counter = counter + 1
           data(i,j) = recv_buffer_2d(counter)/glocal%size
        end do
     end do
     is_get_ok = .true.
  else
     is_get_ok = .false.
  end if
   
end subroutine get_ensemble_mean_2d

  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine bcast_ensemble_int(data)
  implicit none
  integer, intent(INOUT) :: data
  integer :: int_array(1)

  if (glocal%rank == 0) int_array(1) = data

  call MPI_Bcast(int_array, 1, MPI_INT, 0, glocal%comm, ierror)

  data = int_array(1)

end subroutine bcast_ensemble_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine bcast_ensemble_1d(data)
  implicit none
  real(kind=8), intent(INOUT) :: data(:)
  integer :: size1

  size1 = size(data)
  
  call MPI_Bcast(data, size1, MPI_DOUBLE, 0, glocal%comm, ierror)
  
end subroutine bcast_ensemble_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine bcast_ensemble_2d(data)
  implicit none
  real(kind=8), intent(INOUT) :: data(:,:)
  integer :: size1, size2

  size1 = size(data,1)
  size2 = size(data,2)
  
  call MPI_Bcast(data, size1*size2, MPI_DOUBLE, 0, glocal%comm, ierror)
  
end subroutine bcast_ensemble_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_int_1d_local(data,is,ie,dest)
  implicit none
  integer, intent(IN) :: data(:)
  integer, intent(IN) :: is, ie
  integer, intent(IN) :: dest

  integer :: buffer(is:ie)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie) = data(is:ie)
  call MPI_ISEND(buffer,ie-is+1,MPI_INTEGER,dest,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_int_1d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_int_1d_local(data,is,ie,source)
  implicit none
  integer, intent(INOUT) :: data(:)
  integer, intent(IN)  :: is, ie
  integer, intent(IN)  :: source
  
  integer :: buffer(is:ie)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,ie-is+1,MPI_INTEGER,source,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie) = buffer(is:ie)

end subroutine recv_int_1d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_int_2d_local(data,is,ie,js,je,dest_pe)
  implicit none
  integer, intent(IN) :: data(:,:)
  integer, intent(IN) :: is, ie, js, je
  integer, intent(IN) :: dest_pe

  integer :: buffer(is:ie,js:je)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie,js:je) = data(is:ie,js:je)

  call MPI_ISEND(buffer,(ie-is+1)*(je-js+1),MPI_INTEGER,dest_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_int_2d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_int_2d_local(data,is,ie,js,je,source_pe)
  implicit none
  integer, intent(INOUT) :: data(:,:)
  integer, intent(IN) :: is, ie, js, je
  integer, intent(IN) :: source_pe

  integer :: buffer(is:ie,js:je)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,(ie-is+1)*(je-js+1),MPI_INTEGER,source_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie,js:je) = buffer(is:ie,js:je)

end subroutine recv_int_2d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_int_3d_local(data,is,ie,js,je,ks,ke,dest_pe)
  implicit none
  integer, intent(IN) :: data(:,:,:)
  integer, intent(IN) :: is, ie, js, je, ks, ke
  integer, intent(IN) :: dest_pe

  integer :: buffer(is:ie,js:je,ks:ke)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie,js:je,ks:ke) = data(is:ie,js:je,ks:ke)

  call MPI_ISEND(buffer,(ie-is+1)*(je-js+1)*(ke-ks+1),MPI_INTEGER,dest_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_int_3d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_int_3d_local(data,is,ie,js,je,ks,ke,source_pe)
  implicit none
  integer, intent(INOUT) :: data(:,:,:)
  integer, intent(IN) :: is, ie, js, je, ks, ke
  integer, intent(IN) :: source_pe

  integer :: buffer(is:ie,js:je,ks:ke)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,(ie-is+1)*(je-js+1)*(ke-ks+1),MPI_INTEGER,source_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie,js:je,ks:ke) = buffer(is:ie,js:je,ks:ke)

end subroutine recv_int_3d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_real_1d_local(data,is,ie,dest)
  implicit none
  real(kind=4), intent(IN) :: data(:)
  integer, intent(IN) :: is, ie
  integer, intent(IN) :: dest

  real(kind=4) :: buffer(is:ie)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie) = data(is:ie)
  call MPI_ISEND(buffer,ie-is+1,MPI_REAL,dest,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_real_1d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_real_1d_local(data,is,ie,source)
  implicit none
  real(kind=4), intent(INOUT) :: data(:)
  integer, intent(IN)  :: is, ie
  integer, intent(IN)  :: source
  
  real(kind=4) :: buffer(is:ie)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,ie-is+1,MPI_REAL,source,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie) = buffer(is:ie)

end subroutine recv_real_1d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_real_2d_local(data,is,ie,js,je,dest_pe)
  implicit none
  real(kind=4), intent(IN) :: data(:,:)
  integer, intent(IN) :: is, ie, js, je
  integer, intent(IN) :: dest_pe

  real(kind=4) :: buffer(is:ie,js:je)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie,js:je) = data(is:ie,js:je)

  call MPI_ISEND(buffer,(ie-is+1)*(je-js+1),MPI_REAL,dest_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_real_2d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_real_2d_local(data,is,ie,js,je,source_pe)
  implicit none
  real(kind=4), intent(INOUT) :: data(:,:)
  integer, intent(IN) :: is, ie, js, je
  integer, intent(IN) :: source_pe

  real(kind=4) :: buffer(is:ie,js:je)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,(ie-is+1)*(je-js+1),MPI_REAL,source_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie,js:je) = buffer(is:ie,js:je)

end subroutine recv_real_2d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_real_3d_local(data,is,ie,js,je,ks,ke,dest_pe)
  implicit none
  real(kind=4), intent(IN) :: data(:,:,:)
  integer, intent(IN) :: is, ie, js, je, ks, ke
  integer, intent(IN) :: dest_pe

  real(kind=4) :: buffer(is:ie,js:je,ks:ke)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie,js:je,ks:ke) = data(is:ie,js:je,ks:ke)

  call MPI_ISEND(buffer,(ie-is+1)*(je-js+1)*(ke-ks+1),MPI_REAL,dest_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_real_3d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_real_3d_local(data,is,ie,js,je,ks,ke,source_pe)
  implicit none
  real(kind=4), intent(INOUT) :: data(:,:,:)
  integer, intent(IN) :: is, ie, js, je, ks, ke
  integer, intent(IN) :: source_pe

  real(kind=4) :: buffer(is:ie,js:je,ks:ke)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,(ie-is+1)*(je-js+1)*(ke-ks+1),MPI_REAL,source_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie,js:je,ks:ke) = buffer(is:ie,js:je,ks:ke)

end subroutine recv_real_3d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_double_1d_local(data,is,ie,dest)
  implicit none
  real(kind=8), intent(IN) :: data(:)
  integer, intent(IN) :: is, ie
  integer, intent(IN) :: dest

  real(kind=8) :: buffer(is:ie)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie) = data(is:ie)
  call MPI_ISEND(buffer,ie-is+1,MPI_DOUBLE,dest,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_double_1d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_double_1d_local(data,is,ie,source)
  implicit none
  real(kind=8), intent(INOUT) :: data(:)
  integer, intent(IN)  :: is, ie
  integer, intent(IN)  :: source
  
  real(kind=8) :: buffer(is:ie)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,ie-is+1,MPI_DOUBLE,source,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie) = buffer(is:ie)

end subroutine recv_double_1d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_double_2d_local(data,is,ie,js,je,dest_pe)
  implicit none
  real(kind=8), intent(IN) :: data(:,:)
  integer, intent(IN) :: is, ie, js, je
  integer, intent(IN) :: dest_pe

  real(kind=8) :: buffer(is:ie,js:je)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie,js:je) = data(is:ie,js:je)

  call MPI_ISEND(buffer,(ie-is+1)*(je-js+1),MPI_DOUBLE,dest_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_double_2d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_double_2d_local(data,is,ie,js,je,source_pe)
  implicit none
  real(kind=8), intent(INOUT) :: data(:,:)
  integer, intent(IN) :: is, ie, js, je
  integer, intent(IN) :: source_pe

  real(kind=8) :: buffer(is:ie,js:je)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,(ie-is+1)*(je-js+1),MPI_DOUBLE,source_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie,js:je) = buffer(is:ie,js:je)

end subroutine recv_double_2d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_double_3d_local(data,is,ie,js,je,ks,ke,dest_pe)
  implicit none
  real(kind=8), intent(IN) :: data(:,:,:)
  integer, intent(IN) :: is, ie, js, je, ks, ke
  integer, intent(IN) :: dest_pe

  real(kind=8) :: buffer(is:ie,js:je,ks:ke)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  buffer(is:ie,js:je,ks:ke) = data(is:ie,js:je,ks:ke)

  call MPI_ISEND(buffer,(ie-is+1)*(je-js+1)*(ke-ks+1),MPI_DOUBLE,dest_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

end subroutine send_double_3d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_double_3d_local(data,is,ie,js,je,ks,ke,source_pe)
  implicit none
  real(kind=8), intent(INOUT) :: data(:,:,:)
  integer, intent(IN) :: is, ie, js, je, ks, ke
  integer, intent(IN) :: source_pe

  real(kind=8) :: buffer(is:ie,js:je,ks:ke)
  integer :: request
  integer :: status(MPI_STATUS_SIZE)

  call MPI_IRECV(buffer,(ie-is+1)*(je-js+1)*(ke-ks+1),MPI_DOUBLE,source_pe,MPI_MY_TAG,local%comm,request,ierror)
  call MPI_WAIT(request,status,ierror)

  data(is:ie,js:je,ks:ke) = buffer(is:ie,js:je,ks:ke)

end subroutine recv_double_3d_local

!=======+=========+=========+=========+=========+=========+=========+=========+

end module mod_ensemble_base

