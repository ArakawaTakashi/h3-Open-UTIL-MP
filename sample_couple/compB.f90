program compB
  use mpi
  use h3ou_api
  implicit none
  character(len=5)          :: my_name = "compB"
  integer                   :: time_array(6) = [2000, 1, 1, 0, 0, 0]
  integer, parameter        :: delta_t = 3600
  integer                   :: my_comm, my_group, my_size, my_rank
  integer, parameter        :: global_array_size = 20
  integer                   :: array_size   ! local array size
  integer, allocatable      :: my_grid(:)
  integer, allocatable      :: send_index(:)
  integer, allocatable      :: recv_index(:)
  real(kind=8), allocatable :: coef(:)
  integer, allocatable      :: my_array_int(:)
  real(kind=4), allocatable :: my_array_real(:)
  real(kind=8), allocatable :: my_array_double(:)
  integer, allocatable      :: target_array_int(:)
  real(kind=4), allocatable :: target_array_real(:)
  real(kind=8), allocatable :: target_array_double(:)
  logical                   :: is_get_ok
  integer                   :: t
  integer                   :: status(MPI_STATUS_SIZE)
  integer                   :: ierror
  integer                   :: i
  
  ! initialize 
  call h3ou_init(my_name, "coupling.conf")

  ! get mpi parameters
  call h3ou_get_mpi_parameter(my_name, my_comm, my_group, my_size,  my_rank)

  write(0, *) my_name//" my_size = ", my_size, ", my_rank = ", my_rank

  if (mod(global_array_size, my_size) /= 0) then
     write(0, *) "Error, mod(global_array_size, my_size) must be 0"
     stop 9999
  end if
  
  array_size  = int(global_array_size/my_size)

  allocate(my_grid(array_size))
  allocate(my_array_int(array_size))
  allocate(my_array_real(array_size))
  allocate(my_array_double(array_size))
  allocate(target_array_int(array_size))
  allocate(target_array_real(array_size))
  allocate(target_array_double(array_size))

  do i = 1, array_size
     my_grid(i) = array_size * my_rank + i
  end do

  call h3ou_def_grid(my_grid, "compB", "compB_grid1", 20)
  call h3ou_end_grid_def()

  if (my_rank == 0) then
     allocate(send_index(global_array_size))
     allocate(recv_index(global_array_size))
     allocate(coef(global_array_size))
     do i = 1, global_array_size
        send_index(i) = i
        recv_index(i) = i
        coef(i)  = 1.d0
     end do
  else
     allocate(send_index(1))
     allocate(recv_index(1))
     allocate(coef(1))
  end if

  call h3ou_set_interpolation_table("compB", "compA", "compA_grid1", "compB", "compB_grid1", 1)
  
  call h3ou_set_interpolation_table("compB", "compB", "compB_grid1", "compA", "compA_grid1", 1, &
       send_index, recv_index, coef)

  call h3ou_init_time(time_array)
  
  ! set my array
  do i = 1, array_size
     my_array_int(i)    = 20000 + my_rank * 100 + i
     my_array_real(i)   = my_array_int(i)
     my_array_double(i) = my_array_int(i)
  end do

  call h3ou_put_data("compB_var1", my_array_double)

  do t = 1, 24
     call h3ou_set_time(my_name, time_array, delta_t)

     call h3ou_get_data("compA_var1", target_array_double, is_recv_ok = is_get_ok)

     call h3ou_put_data("compB_var1", my_array_double)

     call h3ou_inc_time(time_array, delta_t)
     
  end do

  call h3ou_coupling_end(time_array, .true.)

end program compB
