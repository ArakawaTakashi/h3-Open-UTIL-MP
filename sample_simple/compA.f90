program compA
  use h3ou_api
  implicit none
  
  character(len=5)  :: my_name = "compA"
  integer           :: my_comm, my_group, my_size, my_rank
  integer,parameter :: array_size = 10
  integer           :: my_array_int(array_size)
  real(kind=4)      :: my_array_real(array_size)
  real(kind=8)      :: my_array_double(array_size)
  integer           :: target_array_int(array_size)
  real(kind=4)      :: target_array_real(array_size)
  real(kind=8)      :: target_array_double(array_size)
  integer           :: i
  
  ! initialize 
  call h3ou_init(my_name, "coupling.conf")

  ! get mpi parameters
  call h3ou_get_mpi_parameter(my_name, my_comm, my_group, my_size,  my_rank)

  write(0, *) my_name//" my_size = ", my_size, ", my_rank = ", my_rank

  ! set my array
  do i = 1, array_size
     my_array_int(i)    = 10000 + my_rank * 100 + i
     my_array_real(i)   = my_array_int(i)
     my_array_double(i) = my_array_int(i)
  end do
  
  ! bcast global

  if (my_rank == 0) then
     call h3ou_bcast_global("compA", my_array_int)
  else
     call h3ou_bcast_global("compA", target_array_int)
  end if

  call h3ou_bcast_global("compB", target_array_real)

  call h3ou_bcast_global("compC", target_array_double)
  
  ! bcast model

  if (my_rank == 0) then
     call h3ou_bcast_model("compA", "compB", my_array_int)
  else
     target_array_int = 0
     call h3ou_bcast_model("compA", "compB", target_array_int)
     write(0, *) "compA bcast_model, ", target_array_int
  end if
  
  ! send/recv model

  if (my_rank == 0) then
     call h3ou_send_model("compB", 0, my_array_double)
     call h3ou_recv_model("compB", 0, target_array_real)
     write(0, *) "compA send/recv model, ", target_array_real 
  end if
  
  call h3ou_coupling_end()
end program compA
