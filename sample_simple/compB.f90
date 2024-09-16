program compB
  use h3ou_api
  integer :: my_data = 10
  integer :: target_data
  integer :: recv_array(10)
  
  !call h3ou_init_simple("compB", "LOUD")
  call h3ou_init("compB", "coupling.conf")

  call h3ou_recv("compA", target_data)
  call h3ou_send("compA", my_data)
  write(0, *) "compB target_data = ", target_data

  call h3ou_recv_model_int("compA", 0, recv_array)

  write(0, *) "irecv array = ", recv_array
  
  call h3ou_coupling_end()
  
end program compB
