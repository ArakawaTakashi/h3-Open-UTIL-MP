program compA
  use h3ou_api
  integer :: my_data = 1
  integer :: target_data

  call h3ou_init_simple("compA", "LOUD")
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
  call h3ou_end()
end program compA
