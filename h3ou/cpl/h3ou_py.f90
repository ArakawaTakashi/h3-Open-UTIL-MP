!=======+=========+=========+=========+=========+=========+=========+=========+
!
! h3open-util/mp python api
! these routines are called from python h3open-util/mp api h3opp.py and h3ouc_api.c
!
!=======+=========+=========+=========+=========+=========+=========+=========+
subroutine h3oup_init(comp_name, comp_name_len, config_name, config_name_len) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_init
  implicit none
  character(len=1)    :: comp_name(*)
  integer, intent(IN) :: comp_name_len
  character(len=1)    :: config_name(*)
  integer, intent(IN) :: config_name_len

  call h3ou_init(trim(get_char_str(comp_name, comp_name_len)), &
                 trim(get_char_str(config_name, config_name_len)))

end subroutine h3oup_init

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_get_mpi_parameter(comp_name, comp_name_len, my_comm, my_group, my_size, my_rank) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_get_mpi_parameter
  implicit none
  character(len=1)    :: comp_name(*)
  integer, intent(IN) :: comp_name_len
  integer, intent(INOUT) :: my_comm
  integer, intent(INOUT) :: my_group
  integer, intent(INOUT) :: my_size
  integer, intent(INOUT) :: my_rank

  call h3ou_get_mpi_parameter(trim(get_char_str(comp_name, comp_name_len)), my_comm, my_group, my_size, my_rank)
  
end subroutine h3oup_get_mpi_parameter

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function h3oup_get_my_rank() bind(C)
  use h3ou_api, only : h3ou_get_my_rank
  implicit none

  h3oup_get_my_rank = h3ou_get_my_rank()

end function h3oup_get_my_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function h3oup_get_my_size() bind(C)
  use h3ou_api, only : h3ou_get_my_size
  implicit none

  h3oup_get_my_size = h3ou_get_my_size()

end function h3oup_get_my_size

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_def_grid(grid_index, ngrid, &
                         comp_name, len_comp, &
                         grid_name, len_grid, &
                         nz)  bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_def_grid
  integer, intent(IN) :: ngrid
  integer, intent(IN) :: grid_index(ngrid)
  character(len=1), intent(IN) :: comp_name(*)
  integer, intent(IN) :: len_comp
  character(len=1), intent(IN) :: grid_name(*)
  integer, intent(IN) :: len_grid
  integer, intent(IN) :: nz

 
  call h3ou_def_grid(grid_index, trim(get_char_str(comp_name, len_comp)), &
                    trim(get_char_str(grid_name, len_grid)), nz)

end subroutine h3oup_def_grid

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_end_grid_def() bind(C)
  use h3ou_api, only : h3ou_end_grid_def
  implicit none

  call h3ou_end_grid_def()

end subroutine h3oup_end_grid_def

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_set_interpolation_table(my_comp, mcomp_len, &
                                  send_comp, scomp_len, send_grid, sgrid_len, &
                                  recv_comp, rcomp_len, recv_grid, rgrid_len, &
                                  mapping_tag, &
                                  send_index, recv_index, coef, nindex) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_set_interpolation_table
  implicit none
  character(len=1), intent(IN) :: my_comp(*)            
  integer, intent(IN)          :: mcomp_len
  character(len=1), intent(IN) :: send_comp(*)
  integer, intent(IN)          :: scomp_len
  character(len=1), intent(IN) :: send_grid(*)
  integer, intent(IN)          :: sgrid_len
  character(len=1), intent(IN) :: recv_comp(*)
  integer, intent(IN)          :: rcomp_len
  character(len=1), intent(IN) :: recv_grid(*)
  integer, intent(IN)          :: rgrid_len
  integer, intent(IN)          :: mapping_tag
  integer, intent(IN)          :: nindex
  integer, intent(IN)          :: send_index(nindex)
  integer, intent(IN)          :: recv_index(nindex)
  real(kind=8), intent(IN)     :: coef(nindex)
  
  call h3ou_set_interpolation_table(trim(get_char_str(my_comp, mcomp_len)),   &
                                   trim(get_char_str(send_comp, scomp_len)), &
                                   trim(get_char_str(send_grid, sgrid_len)), &
                                   trim(get_char_str(recv_comp, rcomp_len)), &
                                   trim(get_char_str(recv_grid, rgrid_len)), &
                                   mapping_tag, &
                                   send_index = send_index, &
                                   recv_index = recv_index, &
                                   coef = coef)


end subroutine h3oup_set_interpolation_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_set_interpolation_table_no_index(my_comp, mcomp_len, &
                                  send_comp, scomp_len, send_grid, sgrid_len, &
                                  recv_comp, rcomp_len, recv_grid, rgrid_len, &
                                  mapping_tag) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_set_interpolation_table
  implicit none
  character(len=1), intent(IN) :: my_comp(*)
  integer, intent(IN)          :: mcomp_len
  character(len=1), intent(IN) :: send_comp(*)
  integer, intent(IN)          :: scomp_len
  character(len=1), intent(IN) :: send_grid(*)
  integer, intent(IN)          :: sgrid_len
  character(len=1), intent(IN) :: recv_comp(*)
  integer, intent(IN)          :: rcomp_len
  character(len=1), intent(IN) :: recv_grid(*)
  integer, intent(IN)          :: rgrid_len
  integer, intent(IN)          :: mapping_tag
  
  call h3ou_set_interpolation_table(trim(get_char_str(my_comp, mcomp_len)),   &
                                   trim(get_char_str(send_comp, scomp_len)), &
                                   trim(get_char_str(send_grid, sgrid_len)), &
                                   trim(get_char_str(recv_comp, rcomp_len)), &
                                   trim(get_char_str(recv_grid, rgrid_len)), &
                                   mapping_tag)

end subroutine h3oup_set_interpolation_table_no_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_set_mapping_table(my_comp, mcomp_len, &
                                  send_comp, scomp_len, send_grid, sgrid_len, &
                                  recv_comp, rcomp_len, recv_grid, rgrid_len, &
                                  mapping_tag, &
                                  send_index, recv_index, nindex) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_set_interpolation_table
  implicit none
  character(len=1), intent(IN) :: my_comp(*)            
  integer, intent(IN)          :: mcomp_len
  character(len=1), intent(IN) :: send_comp(*)
  integer, intent(IN)          :: scomp_len
  character(len=1), intent(IN) :: send_grid(*)
  integer, intent(IN)          :: sgrid_len
  character(len=1), intent(IN) :: recv_comp(*)
  integer, intent(IN)          :: rcomp_len
  character(len=1), intent(IN) :: recv_grid(*)
  integer, intent(IN)          :: rgrid_len
  integer, intent(IN)          :: mapping_tag
  integer, intent(IN)          :: nindex
  integer, intent(IN)          :: send_index(nindex)
  integer, intent(IN)          :: recv_index(nindex)
  
  call h3ou_set_interpolation_table(trim(get_char_str(my_comp, mcomp_len)),   &
                                   trim(get_char_str(send_comp, scomp_len)), &
                                   trim(get_char_str(send_grid, sgrid_len)), &
                                   trim(get_char_str(recv_comp, rcomp_len)), &
                                   trim(get_char_str(recv_grid, rgrid_len)), &
                                   mapping_tag, &
                                   send_index = send_index, &
                                   recv_index = recv_index)

end subroutine h3oup_set_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_set_mapping_table_no_index(my_comp, mcomp_len, &
                                  send_comp, scomp_len, send_grid, sgrid_len, &
                                  recv_comp, rcomp_len, recv_grid, rgrid_len, &
                                  mapping_tag) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_set_interpolation_table
  implicit none
  character(len=1), intent(IN) :: my_comp(*)
  integer, intent(IN)          :: mcomp_len
  character(len=1), intent(IN) :: send_comp(*)
  integer, intent(IN)          :: scomp_len
  character(len=1), intent(IN) :: send_grid(*)
  integer, intent(IN)          :: sgrid_len
  character(len=1), intent(IN) :: recv_comp(*)
  integer, intent(IN)          :: rcomp_len
  character(len=1), intent(IN) :: recv_grid(*)
  integer, intent(IN)          :: rgrid_len
  integer, intent(IN)          :: mapping_tag
  
  call h3ou_set_interpolation_table(trim(get_char_str(my_comp, mcomp_len)),   &
                                   trim(get_char_str(send_comp, scomp_len)), &
                                   trim(get_char_str(send_grid, sgrid_len)), &
                                   trim(get_char_str(recv_comp, rcomp_len)), &
                                   trim(get_char_str(recv_grid, rgrid_len)), &
                                   mapping_tag)

end subroutine h3oup_set_mapping_table_no_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_get_num_of_put_data(num_of_pdata) bind(C)
  use h3ou_api, only : h3ou_get_num_of_put_data
  implicit none
  integer, intent(OUT) :: num_of_pdata
  
  num_of_pdata = h3ou_get_num_of_put_data()

end subroutine h3oup_get_num_of_put_data

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_get_put_data_name(pdata_num, data_name, name_len) bind(C)
  use h3ou_api, only : h3ou_get_put_data_name
  implicit none
  integer, intent(IN) :: pdata_num
  character(len=1), intent(INOUT) :: data_name(*)
  integer, intent(INOUT)          :: name_len
  character(len=64) :: dname
  integer :: i
  
  dname = trim(h3ou_get_put_data_name(pdata_num))

  do i = 1, len_trim(dname)
     data_name(i) = dname(i:i)
  end do

  data_name(len_trim(dname)+1) = CHAR(0)

end subroutine h3oup_get_put_data_name

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_get_num_of_get_data(num_of_gdata) bind(C)
  use h3ou_api, only : h3ou_get_num_of_get_data
  implicit none
  integer, intent(OUT) :: num_of_gdata
  
  num_of_gdata = h3ou_get_num_of_get_data()

end subroutine h3oup_get_num_of_get_data

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_get_get_data_name(gdata_num, data_name, name_len) bind(C)
  use h3ou_api, only : h3ou_get_get_data_name
  implicit none
  integer, intent(IN) :: gdata_num
  character(len=1), intent(INOUT) :: data_name(*)
  integer, intent(INOUT)          :: name_len
  character(len=64) :: dname
  integer :: i
  
  dname = trim(h3ou_get_get_data_name(gdata_num))

  do i = 1, len_trim(dname)
     data_name(i) = dname(i:i)
  end do

  data_name(len_trim(dname)+1) = CHAR(0)
  
end subroutine h3oup_get_get_data_name

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_get_vlayer(data_name, name_len, vlayer) bind(C)
  use h3ou_api, only : h3ou_get_vlayer
  use h3ou_py_base, only : get_char_str
  implicit none
  character(len=1), intent(IN) :: data_name(*)
  integer, intent(IN)          :: name_len
  integer, intent(OUT)         :: vlayer

  vlayer = h3ou_get_vlayer(trim(get_char_str(data_name, name_len)))

end subroutine h3oup_get_vlayer


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_init_time(time_array, nsize) bind(C)
  use h3ou_api, only : h3ou_init_time
  implicit none
  integer, intent(IN) :: nsize
  integer, intent(IN) :: time_array(nsize)

  call h3ou_init_time(time_array)

end subroutine h3oup_init_time

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_set_time(my_comp, mcomp_len, time_array, nsize, delta_t) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_set_time
  implicit none
  character(len=1), intent(IN) :: my_comp(*)
  integer, intent(IN)          :: mcomp_len
  integer, intent(IN)          :: nsize
  integer, intent(IN)          :: time_array(nsize)
  integer, intent(IN)          :: delta_t
  
  call h3ou_set_time(trim(get_char_str(my_comp, mcomp_len)), time_array, delta_t)

end subroutine h3oup_set_time

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_put_data_1d(data_name, name_len, data, data_len) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_put_data_1d
  implicit none
  character(len=1), intent(IN) :: data_name(*)
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: data_len
  real(kind=8), intent(INOUT)  :: data(data_len)

  call h3ou_put_data_1d(trim(get_char_str(data_name, name_len)), data)

end subroutine h3oup_put_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_put_data_25d(data_name, name_len, data, data_len1, data_len2) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_put_data_25d
  implicit none
  character(len=1), intent(IN) :: data_name(*)
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: data_len1
  integer, intent(IN)          :: data_len2
  real(kind=8), intent(INOUT)  :: data(data_len1*data_len2)
  real(kind=8) :: data_buffer(data_len2, data_len1)
  integer :: datashape(2)

  datashape(1) = data_len2
  datashape(2) = data_len1
  data_buffer = reshape(data, datashape)
  
  call h3ou_put_data_25d(trim(get_char_str(data_name, name_len)), data_buffer)

end subroutine h3oup_put_data_25d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_get_data_1d(data_name, name_len, data, data_len, is_recv_ok) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_get_data_1d
  implicit none
  character(len=1), intent(IN) :: data_name(*)
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: data_len
  real(kind=8), intent(INOUT)  :: data(data_len)
  integer, intent(INOUT)       :: is_recv_ok
  logical :: get_flag
  
  call h3ou_get_data_1d(trim(get_char_str(data_name, name_len)), data, is_recv_ok = get_flag)

  if (get_flag) then
     is_recv_ok = 1
  else
     is_recv_ok = 0
  end if

end subroutine h3oup_get_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_get_data_25d(data_name, name_len, data, data_len1, data_len2, is_recv_ok) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_get_data_25d
  implicit none
  character(len=1), intent(IN) :: data_name(*)
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: data_len1
  integer, intent(IN)          :: data_len2
  real(kind=8), intent(INOUT)  :: data(data_len1*data_len2)
  integer, intent(INOUT)       :: is_recv_ok
  real(kind=8) :: data_buffer(data_len2, data_len1)
  integer :: datashape(1)
  logical :: get_flag
  
  call h3ou_get_data_25d(trim(get_char_str(data_name, name_len)), data_buffer, is_recv_ok = get_flag)

  datashape(1) = data_len1*data_len2
  data = reshape(data_buffer, datashape)

  if (get_flag) then
     is_recv_ok = 1
  else
     is_recv_ok = 0
  end if

  return
  
end subroutine h3oup_get_data_25d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_coupling_end(time_array, nsize) bind(C)
  use h3ou_api, only : h3ou_coupling_end
  implicit none
  integer, intent(IN) :: nsize
  integer, intent(IN) :: time_array(nsize)

  call h3ou_coupling_end(time_array, .false.)

end subroutine h3oup_coupling_end

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_array_int(my_comp_name, my_name_len, recv_comp_name, recv_name_len, array, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=1), intent(IN) :: my_comp_name
  integer, intent(IN)          :: my_name_len
  character(len=1), intent(IN) :: recv_comp_name
  integer, intent(IN)          :: recv_name_len
  integer, intent(IN)          :: array_size
  integer, intent(IN)          :: array(array_size)

  call h3ou_send_array(trim(get_char_str(my_comp_name, my_name_len)), trim(get_char_str(recv_comp_name, recv_name_len)), array)
  
end subroutine h3oup_send_array_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_array_int(my_comp_name, my_name_len, send_comp_name, send_name_len, array, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=1), intent(IN) :: my_comp_name
  integer, intent(IN)          :: my_name_len
  character(len=1), intent(IN) :: send_comp_name
  integer, intent(IN)          :: send_name_len
  integer, intent(IN)          :: array_size
  integer, intent(INOUT)       :: array(array_size)

  call h3ou_recv_array(trim(get_char_str(my_comp_name, my_name_len)), trim(get_char_str(send_comp_name, send_name_len)), array)
  
end subroutine h3oup_recv_array_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_array_double(my_comp_name, my_name_len, recv_comp_name, recv_name_len, array, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=1), intent(IN) :: my_comp_name
  integer, intent(IN)          :: my_name_len
  character(len=1), intent(IN) :: recv_comp_name
  integer, intent(IN)          :: recv_name_len
  integer, intent(IN)          :: array_size
  real(kind=8), intent(IN)     :: array(array_size)

  call h3ou_send_array(trim(get_char_str(my_comp_name, my_name_len)), trim(get_char_str(recv_comp_name, recv_name_len)), array)
  
end subroutine h3oup_send_array_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_array_double(my_comp_name, my_name_len, send_comp_name, send_name_len, array, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=1), intent(IN) :: my_comp_name
  integer, intent(IN)          :: my_name_len
  character(len=1), intent(IN) :: send_comp_name
  integer, intent(IN)          :: send_name_len
  integer, intent(IN)          :: array_size
  real(kind=8), intent(INOUT)  :: array(array_size)

  call h3ou_recv_array(trim(get_char_str(my_comp_name, my_name_len)), trim(get_char_str(send_comp_name, send_name_len)), array)
  
end subroutine h3oup_recv_array_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_inc_calendar(time_array, nsize, delta_t) bind(C)
  use h3ou_api, only : h3ou_inc_calendar
  integer, intent(IN)    :: nsize
  integer, intent(INOUT) :: time_array(nsize)
  integer, intent(IN)    :: delta_t

  call h3ou_inc_calendar(time_array, delta_t)
  
end subroutine h3oup_inc_calendar

!=======+=========+=========+=========+=========+=========+=========+=========+

!subroutine interpolate_data(recv_model_name, send_model_name, mapping_tag, sn1, sn2, send_data, & 
!                            rn1, rn2, recv_data, num_of_data, tn, exchange_tag)
!  use h3ou_api, only : h3ou_interpolation
!  implicit none
!  character(len=*), intent(IN) :: recv_model_name, send_model_name
!  integer, intent(IN) :: mapping_tag
!  integer, intent(IN) :: sn1, sn2
!  real(kind=8), intent(IN) :: send_data(sn1,sn2)
!  integer, intent(IN) :: rn1, rn2
!  real(kind=8), intent(INOUT) :: recv_data(rn1,rn2)
!  integer, intent(IN) :: num_of_data
!  integer, intent(IN) :: tn
!  integer, intent(IN) :: exchange_tag(tn)

!  call h3ou_interpolation(recv_model_name, send_model_name, mapping_tag, &
!                          send_data, recv_data, num_of_data, exchange_tag(1))

!end subroutine interpolate_data

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+            Simple Mode APIs           +=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_end() bind(C)
  use h3ou_api, only : h3ou_end
  implicit none

  call h3ou_end()

end subroutine h3oup_end

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+            Global Routines            +=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_global_int(source_name, name_len, data, data_len) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_bcast_global
  implicit none
  character(len=1), intent(IN)    :: source_name
  integer, intent(IN)             :: name_len
  integer, intent(IN)             :: data_len
  integer, intent(INOUT)          :: data(data_len)

  call h3ou_bcast_global(trim(get_char_str(source_name, name_len)), data)

end subroutine h3oup_bcast_global_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_global_real(source_name, name_len, data, data_len) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_bcast_global
  implicit none
  character(len=1), intent(IN)    :: source_name
  integer, intent(IN)             :: name_len
  integer, intent(IN)             :: data_len
  real(kind=4), intent(INOUT)     :: data(data_len)

  call h3ou_bcast_global(trim(get_char_str(source_name, name_len)), data)

end subroutine h3oup_bcast_global_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_global_double(source_name, name_len, data, data_len) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_bcast_global
  implicit none
  character(len=1), intent(IN)    :: source_name
  integer, intent(IN)             :: name_len
  integer, intent(IN)             :: data_len
  real(kind=8), intent(INOUT)     :: data(data_len)

  call h3ou_bcast_global(trim(get_char_str(source_name, name_len)), data)

end subroutine h3oup_bcast_global_double

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+           Inter Model Routines        +=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_model_int(source_name, source_len, target_name, target_len, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_bcast_model
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: source_len
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: target_len
  integer, intent(IN)          :: array_size
  integer, intent(INOUT)       :: data(array_size)

  call h3ou_bcast_model(trim(get_char_str(source_name, source_len)), \
                        trim(get_char_str(target_name, target_len)), data)
  
end subroutine h3oup_bcast_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_model_real(source_name, source_len, target_name, target_len, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_bcast_model
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: source_len
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: target_len
  integer, intent(IN)          :: array_size
  real(kind=4), intent(INOUT)  :: data(array_size)

  call h3ou_bcast_model(trim(get_char_str(source_name, source_len)), \
                        trim(get_char_str(target_name, target_len)), data)

                      end subroutine h3oup_bcast_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_model_double(source_name, source_len, target_name, target_len, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_bcast_model
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: source_len
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: target_len
  integer, intent(IN)          :: array_size
  real(kind=8), intent(INOUT)  :: data(array_size)

  call h3ou_bcast_model(trim(get_char_str(source_name, source_len)), \
                        trim(get_char_str(target_name, target_len)), data)
  
end subroutine h3oup_bcast_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_model_int(target_name, name_len, target_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send_model_int
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: target_pe
  integer, intent(IN)          :: array_size
  integer, intent(IN)          :: data(array_size)

  call h3ou_send_model_int(trim(get_char_str(target_name, name_len)), target_pe, data)
  
end subroutine h3oup_send_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_model_real(target_name, name_len, target_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send_model_real
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: target_pe
  integer, intent(IN)          :: array_size
  real(kind=4), intent(IN)     :: data(array_size)

  call h3ou_send_model_real(trim(get_char_str(target_name, name_len)), target_pe, data)
  
end subroutine h3oup_send_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_model_double(target_name, name_len, target_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send_model_double
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: target_pe
  integer, intent(IN)          :: array_size
  real(kind=8), intent(IN)     :: data(array_size)

  call h3ou_send_model_double(trim(get_char_str(target_name, name_len)), target_pe, data)
  
end subroutine h3oup_send_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_model_int(source_name, name_len, source_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv_model_int
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: source_pe
  integer, intent(IN)          :: array_size
  integer, intent(INOUT)       :: data(array_size)

  call h3ou_recv_model_int(trim(get_char_str(source_name, name_len)), source_pe, data)
  
end subroutine h3oup_recv_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_model_real(source_name, name_len, source_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv_model_real
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: source_pe
  integer, intent(IN)          :: array_size
  real(kind=4), intent(INOUT)  :: data(array_size)

  call h3ou_recv_model_real(trim(get_char_str(source_name, name_len)), source_pe, data)
  
end subroutine h3oup_recv_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_model_double(source_name, name_len, source_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv_model_double
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: source_pe
  integer, intent(IN)          :: array_size
  real(kind=8), intent(INOUT)  :: data(array_size)

  call h3ou_recv_model_double(trim(get_char_str(source_name, name_len)), source_pe, data)
  
end subroutine h3oup_recv_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_int_scalar(target_name, name_len, val) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: val
  
  call h3ou_send(trim(get_char_str(target_name, name_len)), val)
  
end subroutine h3oup_send_int_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_real_scalar(target_name, name_len, val) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  real(kind=4), intent(IN)     :: val
  
  call h3ou_send(trim(get_char_str(target_name, name_len)), val)
  
end subroutine h3oup_send_real_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_double_scalar(target_name, name_len, val) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  real(kind=8), intent(IN)     :: val
  
  call h3ou_send(trim(get_char_str(target_name, name_len)), val)
  
end subroutine h3oup_send_double_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_int_array(target_name, name_len, val, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: array_size
  integer, intent(IN)          :: val(array_size)
  
  call h3ou_send_array(trim(get_char_str(target_name, name_len)), val)
  
end subroutine h3oup_send_int_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_real_array(target_name, name_len, val, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: array_size
  real(kind=4), intent(IN)     :: val(array_size)
  
  call h3ou_send_array(trim(get_char_str(target_name, name_len)), val)
  
end subroutine h3oup_send_real_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_double_array(target_name, name_len, val, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: array_size
  real(kind=8), intent(IN)     :: val(array_size)
  
  call h3ou_send_array(trim(get_char_str(target_name, name_len)), val)
  
end subroutine h3oup_send_double_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_int_scalar(source_name, name_len, val) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(INOUT)       :: val


  call h3ou_recv(trim(get_char_str(source_name, name_len)), val)
  
end subroutine h3oup_recv_int_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_real_scalar(source_name, name_len, val) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  real(kind=4), intent(INOUT)  :: val
  
  call h3ou_recv(trim(get_char_str(source_name, name_len)), val)
  
end subroutine h3oup_recv_real_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_double_scalar(source_name, name_len, val) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  real(kind=8), intent(INOUT)  :: val
  
  call h3ou_recv(trim(get_char_str(source_name, name_len)), val)
  
end subroutine h3oup_recv_double_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_int_array(source_name, name_len, val, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: array_size
  integer, intent(INOUT)       :: val(array_size)
  
  call h3ou_recv_array(trim(get_char_str(source_name, name_len)), val)
  
end subroutine h3oup_recv_int_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_real_array(source_name, name_len, val, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: array_size
  real(kind=4), intent(INOUT)  :: val(array_size)
  
  call h3ou_recv_array(trim(get_char_str(source_name, name_len)), val)
  
end subroutine h3oup_recv_real_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_double_array(source_name, name_len, val, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: array_size
  real(kind=8), intent(INOUT)  :: val(array_size)
  
  call h3ou_recv_array(trim(get_char_str(source_name, name_len)), val)
  
end subroutine h3oup_recv_double_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_isend_model_int(target_name, name_len, target_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_isend_model_int
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: target_pe
  integer, intent(IN)          :: array_size
  integer, intent(IN)          :: data(array_size)

  call h3ou_isend_model_int(trim(get_char_str(target_name, name_len)), target_pe, data)
  
end subroutine h3oup_isend_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_isend_model_real(target_name, name_len, target_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_isend_model_real
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: target_pe
  integer, intent(IN)          :: array_size
  real(kind=4), intent(IN)     :: data(array_size)

  call h3ou_isend_model_real(trim(get_char_str(target_name, name_len)), target_pe, data)
  
end subroutine h3oup_isend_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_isend_model_double(target_name, name_len, target_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_isend_model_double
  implicit none
  character(len=1), intent(IN) :: target_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: target_pe
  integer, intent(IN)          :: array_size
  real(kind=8), intent(IN)     :: data(array_size)

  call h3ou_isend_model_double(trim(get_char_str(target_name, name_len)), target_pe, data)
  
end subroutine h3oup_isend_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_irecv_model_int(source_name, name_len, source_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_irecv_model_int
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: source_pe
  integer, intent(IN)          :: array_size
  integer, intent(INOUT)       :: data(array_size)

  call h3ou_irecv_model_int(trim(get_char_str(source_name, name_len)), source_pe, data)
  
end subroutine h3oup_irecv_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_irecv_model_real(source_name, name_len, source_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_irecv_model_real
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: source_pe
  integer, intent(IN)          :: array_size
  real(kind=4), intent(INOUT)  :: data(array_size)

  call h3ou_irecv_model_real(trim(get_char_str(source_name, name_len)), source_pe, data)
  
end subroutine h3oup_irecv_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_irecv_model_double(source_name, name_len, source_pe, data, array_size) bind(C)
  use h3ou_py_base, only : get_char_str
  use h3ou_api, only : h3ou_irecv_model_double
  implicit none
  character(len=1), intent(IN) :: source_name
  integer, intent(IN)          :: name_len
  integer, intent(IN)          :: source_pe
  integer, intent(IN)          :: array_size
  real(kind=8), intent(INOUT)  :: data(array_size)

  call h3ou_irecv_model_double(trim(get_char_str(source_name, name_len)), source_pe, data)
  
end subroutine h3oup_irecv_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_isend_waitall() bind(C)
  use h3ou_api, only : h3ou_isend_waitall
  implicit none

  call h3ou_isend_waitall()
  
end subroutine h3oup_isend_waitall

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_irecv_waitall() bind(C)
  use h3ou_api, only : h3ou_irecv_waitall
  implicit none

  call h3ou_irecv_waitall()
  
end subroutine h3oup_irecv_waitall

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_local_int(source_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_bcast_local
  implicit none
  integer, intent(IN)    :: source_pe
  integer, intent(IN)    :: array_size
  integer, intent(INOUT) :: data(array_size)

  call h3ou_bcast_local(source_pe, data)

end subroutine h3oup_bcast_local_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_local_real(source_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_bcast_local
  implicit none
  integer, intent(IN)         :: source_pe
  integer, intent(IN)         :: array_size
  real(kind=4), intent(INOUT) :: data(array_size)

  call h3ou_bcast_local(source_pe, data)

end subroutine h3oup_bcast_local_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_bcast_local_double(source_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_bcast_local
  implicit none
  integer, intent(IN)         :: source_pe
  integer, intent(IN)         :: array_size
  real(kind=8), intent(INOUT) :: data(array_size)

  call h3ou_bcast_local(source_pe, data)

end subroutine h3oup_bcast_local_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_local_int(target_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_send_local_all
  integer, intent(IN) :: target_pe
  integer, intent(IN) :: array_size
  integer, intent(IN) :: data(array_size)

  call h3ou_send_local_all(target_pe, data)

end subroutine h3oup_send_local_int
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_local_real(target_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_send_local_all
  integer, intent(IN) :: target_pe
  integer, intent(IN) :: array_size
  real(kind=4), intent(IN) :: data(array_size)

  call h3ou_send_local_all(target_pe, data)

end subroutine h3oup_send_local_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_send_local_double(target_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_send_local_all
  integer, intent(IN) :: target_pe
  integer, intent(IN) :: array_size
  real(kind=8), intent(IN) :: data(array_size)

  call h3ou_send_local_all(target_pe, data)

end subroutine h3oup_send_local_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_local_int(source_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_recv_local_all
  integer, intent(IN) :: source_pe
  integer, intent(IN) :: array_size
  integer, intent(INOUT) :: data(array_size)

  call h3ou_recv_local_all(source_pe, data)

end subroutine h3oup_recv_local_int
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_local_real(source_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_recv_local_all
  integer, intent(IN) :: source_pe
  integer, intent(IN) :: array_size
  real(kind=4), intent(INOUT) :: data(array_size)

  call h3ou_recv_local_all(source_pe, data)

end subroutine h3oup_recv_local_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3oup_recv_local_double(source_pe, data, array_size) bind(C)
  use h3ou_api, only : h3ou_recv_local_all
  integer, intent(IN) :: source_pe
  integer, intent(IN) :: array_size
  real(kind=8), intent(INOUT) :: data(array_size)

  call h3ou_recv_local_all(source_pe, data)

end subroutine h3oup_recv_local_double
  
!=======+=========+=========+=========+=========+=========+=========+=========+
