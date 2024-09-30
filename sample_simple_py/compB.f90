program compB
  use h3ou_api
  implicit none
  character(len=5)  :: my_name = "compB"
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
     my_array_int(i)    = 20000 + my_rank * 100 + i
     my_array_real(i)   = my_array_int(i)
     my_array_double(i) = my_array_int(i)
  end do

  ! bcast global

  call h3ou_bcast_global("compA", target_array_int)

  if (my_rank == 0) then
     call h3ou_bcast_global("compB", my_array_real)
  else
     call h3ou_bcast_global("compB", target_array_real)
  end if
  
  call h3ou_bcast_global("compD", target_array_double)

  ! bcast model

  target_array_int = 0
  call h3ou_bcast_model("compA", "compB", target_array_int)
  write(0, *) "compB bcast model, ", target_array_int

  call h3ou_bcast_model("compB", "compD", my_array_real)
  call h3ou_bcast_model("compD", "compB", target_array_double)
  
  ! send/recv model

  if (my_rank == 0) then
     call h3ou_recv_model("compA", 0, target_array_double)
     call h3ou_send_model("compA", 0, my_array_real)
     write(0, *) "compB send/recv model, ", target_array_double
     call h3ou_send_model("compD", 0, my_array_real)
     call h3ou_recv_model("compD", 0, target_array_double)
  end if
  

  call h3ou_coupling_end()
  
end program compB
