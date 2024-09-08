!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
module jcup_grid_send ! send side interpolation 
  private

!--------------------------------   public  ----------------------------------!

public :: cal_multi_rank_intpl_grid        ! subroutine (send_pe, recv_pe, send_grid, recv_grid)
public :: send_intpl_grid_info             ! subroutine (comp_id)
public :: recv_intpl_grid_info             ! subroutine (comp_id)

!--------------------------------   private  ---------------------------------!

type multi_rank_type
  integer :: recv_grid_index     
  integer :: start_index, end_index        ! global operation index 
  integer :: num_of_send_rank
  integer, pointer :: send_rank(:)
  integer, pointer :: send_grid(:)
  logical :: is_my_intpl = .false.
  integer :: intpl_pe 
end type multi_rank_type

integer :: num_of_multi_rank_intpl = 0 
type(multi_rank_type), pointer :: multi_rank_intpl(:)

integer :: num_of_my_intpl = 0
type(multi_rank_type), pointer :: my_intpl(:)

contains


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine cal_multi_rank_intpl_grid(send_pe, recv_pe, send_grid, recv_grid)
  implicit none
  integer, intent(IN) :: send_pe(:)    ! rank of send component
  integer, intent(IN) :: recv_pe(:)    ! rank of recv component
  integer, intent(IN) :: send_grid(:)  ! send grid index
  integer, intent(IN) :: recv_grid(:)  ! recv grid index
  integer :: recv_index, start_index, end_index
  integer :: num_of_operation
  integer :: n_intpl
  integer :: i
  
  num_of_operation = size(recv_grid)
  recv_index = recv_grid(1)
  start_index = 1
  end_index   = 0
  n_intpl = 0  
  do i = 2, num_of_operation
     if (recv_grid(i) > recv_index) then ! new recv index
        end_index = i-1
        if (minval(send_pe(start_index:end_index)) /= maxval(send_pe(start_index:end_index))) then
          n_intpl = n_intpl + 1
        end if
        start_index = i
        recv_index = recv_grid(i)
     end if
  end do

  num_of_multi_rank_intpl = n_intpl

  allocate(multi_rank_intpl(num_of_multi_rank_intpl))

  recv_index = recv_grid(1)
  start_index = 1
  end_index   = 0
  n_intpl = 0  
  do i = 2, num_of_operation
     if (recv_grid(i) > recv_index) then ! new recv index
        end_index = i-1
        if (minval(send_pe(start_index:end_index)) /= maxval(send_pe(start_index:end_index))) then
           n_intpl = n_intpl + 1
           multi_rank_intpl(n_intpl)%recv_grid_index = recv_index
           multi_rank_intpl(n_intpl)%num_of_send_rank = end_index - start_index + 1
           multi_rank_intpl(n_intpl)%start_index = start_index
           multi_rank_intpl(n_intpl)%end_index   = end_index
           allocate(multi_rank_intpl(n_intpl)%send_rank(end_index - start_index + 1))
           allocate(multi_rank_intpl(n_intpl)%send_grid(end_index - start_index + 1))
           multi_rank_intpl(n_intpl)%send_rank(:) = send_pe(start_index:end_index)           
           multi_rank_intpl(n_intpl)%send_grid(:) = send_grid(start_index:end_index)           
        end if
        start_index = i
        recv_index = recv_grid(i)
     end if
  end do

  
end subroutine cal_multi_rank_intpl_grid

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_intpl_grid_info(comp_id)
  use jcup_mpi_lib, only : jml_BcastLocal
  implicit none
  integer, intent(IN) :: comp_id
  integer :: send_counter(1) ! size of send data
  integer, pointer :: send_buffer(:)
  integer :: counter
  integer :: i, j
  
  
  send_counter(1) = 1
    
  do j = 1, num_of_multi_rank_intpl
    send_counter(1) = send_counter(1) + multi_rank_intpl(j)%num_of_send_rank * 2 + 4
  end do

  call jml_BcastLocal(comp_id, send_counter, 1, 1)

  allocate(send_buffer(send_counter(1)))

  send_buffer(1) = num_of_multi_rank_intpl
  counter = 1

  do j = 1, num_of_multi_rank_intpl
     counter = counter + 1
     send_buffer(counter) = multi_rank_intpl(j)%recv_grid_index
     counter = counter + 1
     send_buffer(counter) = multi_rank_intpl(j)%start_index
     counter = counter + 1
     send_buffer(counter) = multi_rank_intpl(j)%end_index
     counter = counter + 1
     send_buffer(counter) = multi_rank_intpl(j)%num_of_send_rank

     do i = 1, multi_rank_intpl(j)%num_of_send_rank
        counter = counter + 1
        send_buffer(counter) = multi_rank_intpl(j)%send_rank(i)
     end do

     do i = 1, multi_rank_intpl(j)%num_of_send_rank
        counter = counter + 1
        send_buffer(counter) = multi_rank_intpl(j)%send_grid(i)
     end do

  end do

  call jml_BcastLocal(comp_id, send_buffer, 1, send_counter(1))

  deallocate(send_buffer)
  
  call cal_my_intpl(comp_id)
  
end subroutine send_intpl_grid_info

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_intpl_grid_info(comp_id)
  use jcup_mpi_lib, only : jml_BcastLocal
  implicit none
  integer, intent(IN) :: comp_id
  integer :: recv_counter(1)
  integer, pointer :: recv_buffer(:)
  integer :: counter
  integer :: i, j
  
  call jml_BcastLocal(comp_id, recv_counter, 1, 1)

  allocate(recv_buffer(recv_counter(1)))

  call jml_BcastLocal(comp_id, recv_buffer, 1, recv_counter(1))

  num_of_multi_rank_intpl = recv_buffer(1)

  allocate(multi_rank_intpl(num_of_multi_rank_intpl))
  
  counter = 1

  do j = 1, num_of_multi_rank_intpl
     counter = counter + 1
     multi_rank_intpl(j)%recv_grid_index  = recv_buffer(counter)
     counter = counter + 1
     multi_rank_intpl(j)%start_index      = recv_buffer(counter)
     counter = counter + 1
     multi_rank_intpl(j)%end_index        = recv_buffer(counter)
     counter = counter + 1
     multi_rank_intpl(j)%num_of_send_rank = recv_buffer(counter)

     allocate(multi_rank_intpl(j)%send_rank(multi_rank_intpl(j)%num_of_send_rank))
     allocate(multi_rank_intpl(j)%send_grid(multi_rank_intpl(j)%num_of_send_rank))
     
     do i = 1, multi_rank_intpl(j)%num_of_send_rank
        counter = counter + 1
        multi_rank_intpl(j)%send_rank(i)  = recv_buffer(counter)
     end do

     do i = 1, multi_rank_intpl(j)%num_of_send_rank
        counter = counter + 1
        multi_rank_intpl(j)%send_grid(i)  = recv_buffer(counter)
     end do

  end do

  deallocate(recv_buffer)
  
  call cal_my_intpl(comp_id)
  
end subroutine recv_intpl_grid_info

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine cal_my_intpl(comp_id)
  use jcup_mpi_lib, only : jml_GetMyrank
  implicit none
  integer, intent(IN) :: comp_id
  logical, pointer :: my_intpl_flag(:)
  integer :: my_pe ! my rank + 1
  integer :: i, j
 
  my_pe = jml_GetMyrank(comp_id) + 1

  allocate(my_intpl_flag(num_of_multi_rank_intpl))
  my_intpl_flag(:) = .false.
  num_of_my_intpl = 0
  
  do j = 1, num_of_multi_rank_intpl
     do i = 1, multi_rank_intpl(j)%num_of_send_rank
        if (my_pe == multi_rank_intpl(j)%send_rank(i)) then
           num_of_my_intpl = num_of_my_intpl + 1
           my_intpl_flag(num_of_my_intpl) = .true.
           exit
        end if
     end do
  end do


  allocate(my_intpl(num_of_my_intpl))
  
  num_of_my_intpl = 0

  do j = 1, num_of_multi_rank_intpl
     if (my_intpl_flag(j)) then
        num_of_my_intpl = num_of_my_intpl + 1
        my_intpl(num_of_my_intpl)%recv_grid_index  = multi_rank_intpl(j)%recv_grid_index
        my_intpl(num_of_my_intpl)%start_index      = multi_rank_intpl(j)%start_index
        my_intpl(num_of_my_intpl)%end_index        = multi_rank_intpl(j)%end_index
        my_intpl(num_of_my_intpl)%num_of_send_rank = multi_rank_intpl(j)%num_of_send_rank
        allocate(my_intpl(num_of_my_intpl)%send_rank(my_intpl(num_of_my_intpl)%num_of_send_rank))
        allocate(my_intpl(num_of_my_intpl)%send_grid(my_intpl(num_of_my_intpl)%num_of_send_rank))
        my_intpl(num_of_my_intpl)%send_rank(:) = multi_rank_intpl(j)%send_rank(:)
        my_intpl(num_of_my_intpl)%send_grid(:) = multi_rank_intpl(j)%send_grid(:)
     end if
  end do

  do j = 1, num_of_multi_rank_intpl
     deallocate(multi_rank_intpl(j)%send_rank)
     deallocate(multi_rank_intpl(j)%send_grid)
  end do

  deallocate(multi_rank_intpl)
     
end subroutine cal_my_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_my_intpl_flag(my_intpl, my_pe)
  use jcup_utils, only : sort_int_1d
  implicit none
  type(multi_rank_type), intent(INOUT) :: my_intpl
  integer, intent(IN) :: my_pe
  integer, pointer :: sort_array(:)
  integer, pointer :: pe_hist(:), pe_num(:)
  integer :: max_num, max_pe
  integer :: counter
  integer :: i
  
  allocate(sort_array(my_intpl%num_of_send_rank))
  sort_array(:) = my_intpl%send_rank(:)

  call sort_int_1d(my_intpl%num_of_send_rank, sort_array)

  counter = 1
  do i = 2, my_intpl%num_of_send_rank
     if (sort_array(i) > sort_array(i-1)) counter = counter + 1
  end do
  
  allocate(pe_hist(counter))
  allocate(pe_num(counter))
  pe_hist(:) = 0
  counter = 1
  pe_hist(1) = 1
  pe_num(1)  = sort_array(1)
  do i = 2, my_intpl%num_of_send_rank
     if (sort_array(i) > sort_array(i-1)) then
        counter = counter + 1
        pe_num(counter) = sort_array(i)
     end if
     pe_hist(counter) = pe_hist(counter) + 1
  end do

  max_num = 0

  do i = size(pe_hist), 1, -1
     if (max_num <= pe_hist(i)) then
        max_num = pe_hist(i)
        max_pe  = pe_num(i)
     end if
  end do

  if (max_pe == my_pe) my_intpl%is_my_intpl = .true.
  my_intpl%intpl_pe = max_pe
  
  deallocate(pe_hist, pe_num)
  deallocate(sort_array)
  
end subroutine set_my_intpl_flag

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jcup_grid_send
