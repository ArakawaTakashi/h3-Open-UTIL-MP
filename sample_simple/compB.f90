program compB
  use h3ou_api
  integer :: my_data = 10
  integer :: target_data
  
  call h3ou_init_simple("compB", "LOUD")

  call h3ou_recv("compA", target_data)
  call h3ou_send("compA", my_data)
  write(0, *) "compB target_data = ", target_data

  call h3ou_end()
  
end program compB
