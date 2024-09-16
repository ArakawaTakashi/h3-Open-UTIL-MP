program compA
  use h3ou_api
  integer :: my_data = 1
  integer :: target_data
  integer :: my_array(10)
  
  !call h3ou_init_simple("compA", "LOUD")
  call h3ou_init("compA", "coupling.conf")
  my_data = 1
  call h3ou_send("compB", my_data)
  call h3ou_recv("compB", target_data)
  my_data = 2
  call h3ou_send("compC", my_data)
  call h3ou_recv("compC", target_data)
  my_data = 3
  call h3ou_send("compD", my_data)
  call h3ou_recv("compD", target_data)
  write(0, *) "compA target_data = ", target_data

  my_array(:) = 20
  
  call h3ou_send_model_int("compB", 0, my_array)
  call h3ou_send_model_int("compD", 1, my_array)
  
  call h3ou_coupling_end()
end program compA
