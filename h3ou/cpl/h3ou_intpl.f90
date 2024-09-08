!=======+=========+=========+=========+=========+=========+=========+=========+
!> standard interpolation moodule
module h3ou_intpl
  implicit none
  private
!--------------------------------   public  ----------------------------------!

  public :: init_interpolation ! subroutine ()
  public :: set_interpolation ! subroutine (send_comp, send_grid, recv_comp, recv_grid, mapping_tag, coef)
  public :: interpolation ! subroutine (recv_comp, send_comp, mapping_tag, send_data, recv_data, num_of_data, exchange_tag)

  public :: h3ou_intpl_type
  public :: make_parallel_intpl
  
!--------------------------------  private  ----------------------------------!

integer, parameter :: MAX_CHAR_LEN = 14

type intpl_parallel_type
   integer :: num_of_operation
   integer, pointer :: send_data_index(:)
   integer, pointer :: recv_data_index(:)
   real(kind=8), pointer :: coef(:)
end type intpl_parallel_type

type h3ou_intpl_type
  character(len=MAX_CHAR_LEN) :: send_comp
  character(len=MAX_CHAR_LEN) :: send_grid
  character(len=MAX_CHAR_LEN) :: recv_comp ! my_component
  character(len=MAX_CHAR_LEN) :: recv_grid ! my_grid 
  integer :: mapping_tag
  integer :: num_of_operation
  integer, pointer :: operation_index(:)
  integer, pointer :: send_data_index(:)
  integer, pointer :: recv_data_index(:)
  integer, pointer :: send_coef_index(:)
  integer, pointer :: recv_coef_index(:)
  integer, pointer :: global_index_of_local_send_coef(:)
  integer, pointer :: global_index_of_local_recv_coef(:)
  real(kind=8), pointer :: coef(:) => null() ! size(num_of_operation)
  real(kind=8), pointer :: coef_vec(:,:) => null() ! size(num_of_operation, 2)

  integer :: num_of_parallel_intpl
  type(intpl_parallel_type), pointer :: parallel_intpl(:)
  
  type(h3ou_intpl_type), pointer :: next_ptr
end type

integer :: num_of_interpolation = 0
type(h3ou_intpl_type), pointer :: start_ptr => null()
type(h3ou_intpl_type), pointer :: intpl_ptr => null()

contains

!=======+=========+=========+=========+=========+=========+=========+=========+
!>
!! initialize interpolation module
!! \protected
subroutine init_interpolation()
  implicit none

  num_of_interpolation = 0
  start_ptr => null()
  intpl_ptr => null()

end subroutine init_interpolation

!=======+=========+=========+=========+=========+=========+=========+=========+
!>
!! set interpolation information
!! \protected
subroutine set_interpolation(send_comp, send_grid, recv_comp, recv_grid, mapping_tag,  coef, coef_cos, coef_sin)
  implicit none
  character(len=*), intent(IN) :: send_comp   !< name of send component
  character(len=*), intent(IN) :: send_grid   !< name of send grid
  character(len=*), intent(IN) :: recv_comp   !< name of recv component
  character(len=*), intent(IN) :: recv_grid   !< name of recv grid
  integer, intent(IN)          :: mapping_tag !< tag of remapping table
  real(kind=8), optional, intent(IN) :: coef(:) !< interpolation coefficient coef(num_of_global_operation)
  real(kind=8), optional, intent(IN) :: coef_cos(:), coef_sin(:)
  type(h3ou_intpl_type), pointer :: before_ptr

  if (.not.associated(start_ptr)) then
    allocate(start_ptr)
    intpl_ptr => start_ptr
  else

    intpl_ptr => start_ptr
    before_ptr => null()

    do while(associated(intpl_ptr))
       before_ptr => intpl_ptr
       intpl_ptr => intpl_ptr%next_ptr
    end do
  
    allocate(intpl_ptr)
    before_ptr%next_ptr => intpl_ptr
  end if
  
  call set_intpl(intpl_ptr, send_comp, send_grid, recv_comp, recv_grid, mapping_tag, coef, coef_cos, coef_sin)

  call make_parallel_intpl(intpl_ptr)
  
  intpl_ptr%next_ptr => null()

  num_of_interpolation = num_of_interpolation + 1
  
end subroutine set_interpolation

!=======+=========+=========+=========+=========+=========+=========+=========+
!> set intarpolation
!! \private
subroutine set_intpl(current_ptr, send_comp, send_grid, recv_comp, recv_grid, mapping_tag, coef, coef_cos, coef_sin)
  use jcup_interface, only : jcup_get_local_operation_index, jcup_set_local_coef, OPERATION_COEF
  implicit none
  type(h3ou_intpl_type), pointer :: current_ptr
  character(len=*), intent(IN) :: send_comp
  character(len=*), intent(IN) :: send_grid
  character(len=*), intent(IN) :: recv_comp
  character(len=*), intent(IN) :: recv_grid
  integer, intent(IN) :: mapping_tag
  real(kind=8), optional, intent(IN) :: coef(:) ! (num_of_coef, num_of_global_iteration)
  real(kind=8), optional, intent(IN) :: coef_cos(:), coef_sin(:)
  character(len=256) :: log_str

  current_ptr%send_comp = send_comp
  current_ptr%send_grid = send_grid
  current_ptr%recv_comp = recv_comp
  current_ptr%recv_grid = recv_grid
  current_ptr%mapping_tag = mapping_tag
  current_ptr%next_ptr => null()

  call jcup_get_local_operation_index(recv_comp, send_comp, mapping_tag, &
                                      current_ptr%num_of_operation, current_ptr%operation_index, &
                                      current_ptr%send_data_index, current_ptr%recv_data_index, &
                                      current_ptr%send_coef_index, current_ptr%recv_coef_index) 

  current_ptr%send_comp = send_comp
  current_ptr%send_grid = send_grid
  current_ptr%recv_comp = recv_comp
  current_ptr%recv_grid = recv_grid
  current_ptr%mapping_tag = mapping_tag

  if (present(coef)) then
    allocate(current_ptr%coef(current_ptr%num_of_operation))

    call jcup_set_local_coef(recv_comp, send_comp, mapping_tag, coef, current_ptr%coef, OPERATION_COEF)

  else
    !allocate(current_ptr%coef(current_ptr%num_of_operation, 3))
    !current_ptr%coef(:,:) = 1.d0
     current_ptr%coef => null()
  end if

  if (present(coef_cos)) then

    allocate(current_ptr%coef_vec(current_ptr%num_of_operation, 2))
    call jcup_set_local_coef(recv_comp, send_comp, mapping_tag, coef_cos, current_ptr%coef_vec(:,1), OPERATION_COEF)
    call jcup_set_local_coef(recv_comp, send_comp, mapping_tag, coef_sin, current_ptr%coef_vec(:,2), OPERATION_COEF)

  else
     current_ptr%coef_vec => null()
  end if

end subroutine set_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_parallel_intpl(current_ptr)
  use jcup_utils, only : sort_int_1d
  implicit none
  type(h3ou_intpl_type), pointer :: current_ptr
  integer, pointer :: recv_index(:)
  integer, pointer :: conv_table(:)
  integer, pointer :: parallel_index(:,:)
  integer, pointer :: parallel_table(:,:)
  integer, pointer :: index_pointer(:)
  integer :: current_index
  integer :: index_counter
  integer :: max_counter
  integer :: diff_counter
  logical :: is_coef_exist
  integer :: i, j

  allocate(recv_index(current_ptr%num_of_operation))
  allocate(conv_table(current_ptr%num_of_operation))

  do i = 1, current_ptr%num_of_operation
     recv_index(i) = current_ptr%recv_data_index(i)
     conv_table(i) = i
  end do
  
  call  sort_int_1d(current_ptr%num_of_operation, recv_index, conv_table)

  current_index = recv_index(1)
  index_counter = 1
  max_counter = 1
  diff_counter = 1
  do i = 2, current_ptr%num_of_operation -1
     if (current_index == recv_index(i)) then
        index_counter = index_counter + 1
     else
        max_counter = max(max_counter, index_counter)
        index_counter = 1
        diff_counter = diff_counter + 1
     end if
     current_index = recv_index(i)
  end do

  if (current_index == recv_index(current_ptr%num_of_operation)) then
     index_counter = index_counter + 1
     max_counter = max(max_counter, index_counter)
  end if

  if (current_index /= recv_index(current_ptr%num_of_operation)) then
     diff_counter = diff_counter + 1
  end if
  

  allocate(parallel_index(diff_counter, max_counter))
  allocate(parallel_table(diff_counter, max_counter))
  allocate(index_pointer(max_counter))
  
  parallel_index(:,:) = -255
  index_pointer(:) = 1

  current_index = recv_index(1)
  max_counter = 1
  parallel_index(1, 1) = recv_index(1)
  parallel_table(1, 1) = conv_table(1)
  do i = 2, current_ptr%num_of_operation
     if (current_index == recv_index(i)) then ! same number
        max_counter = max_counter + 1
        parallel_index(index_pointer(max_counter), max_counter) = recv_index(i)
        parallel_table(index_pointer(max_counter), max_counter) = conv_table(i)
        index_pointer(max_counter) = index_pointer(max_counter) + 1
     else
        max_counter = 1
        index_pointer(max_counter) = index_pointer(max_counter) + 1
        parallel_index(index_pointer(max_counter), max_counter) = recv_index(i)
        parallel_table(index_pointer(max_counter), max_counter) = conv_table(i)
     end if
     current_index = recv_index(i)
  end do

  current_ptr%num_of_parallel_intpl = size(parallel_index, 2)
  allocate(current_ptr%parallel_intpl(current_ptr%num_of_parallel_intpl))

  is_coef_exist = associated(current_ptr%coef)
  
  do j = 1, current_ptr%num_of_parallel_intpl

     index_counter = 0
     do i = 1, size(parallel_index, 1)
        if (parallel_index(i, j) > 0) then
           index_counter = index_counter + 1
        end if
     end do
     
     current_ptr%parallel_intpl(j)%num_of_operation = index_counter
     allocate(current_ptr%parallel_intpl(j)%send_data_index(index_counter))
     allocate(current_ptr%parallel_intpl(j)%recv_data_index(index_counter))
     if (is_coef_exist) allocate(current_ptr%parallel_intpl(j)%coef(index_counter))

     index_counter = 0
     do i = 1, size(parallel_index, 1)
        if (parallel_index(i, j) > 0) then
           index_counter = index_counter + 1
           current_ptr%parallel_intpl(j)%recv_data_index(index_counter) = parallel_index(i, j)
           current_ptr%parallel_intpl(j)%send_data_index(index_counter) = current_ptr%send_data_index(parallel_table(i,j))
           if (is_coef_exist) current_ptr%parallel_intpl(j)%coef(index_counter) = current_ptr%coef(parallel_table(i,j))
        end if
     end do

  end do

  deallocate(recv_index, conv_table)
  deallocate(parallel_index, parallel_table)
  deallocate(index_pointer)

  !do j = 1, current_ptr%num_of_parallel_intpl
  !   write(0, *) current_ptr%parallel_intpl(j)%recv_data_index(:)
  !end do
  
  !do j = 1, current_ptr%num_of_parallel_intpl
  !   write(0, *) current_ptr%parallel_intpl(j)%send_data_index(:)
  !end do

  !do j = 1, current_ptr%num_of_parallel_intpl
  !   if (is_coef_exist) write(0, *) current_ptr%parallel_intpl(j)%coef(:)
  !end do
  
  return
  
end subroutine make_parallel_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+
!> interpolate 3D vector value
!! \private
subroutine interpolate_vector(send_data, recv_data, num_of_data, &
                                 num_of_operation, send_data_index, recv_data_index, coef1, coef2)
  !use jcup_interface, only : jcup_error
  implicit none
  real(kind=8), intent(IN)    :: send_data(:,:)
  real(kind=8), intent(INOUT) :: recv_data(:,:)
  integer, intent(IN)      :: num_of_data
  integer, intent(IN)      :: num_of_operation
  integer, intent(IN)      :: send_data_index(:)
  integer, intent(IN)      :: recv_data_index(:)
  real(Kind=8), optional, intent(IN) :: coef1(:)
  real(Kind=8), optional, intent(IN) :: coef2(:)
  integer :: i, send_grid, recv_grid
  integer :: k, num_of_vgrid

  !write(6, *) "moj_interpolate_3d_vector IN ", size(send_data, 1), size(send_data, 2), 
  ! num_of_data, num_of_operation, send_data(:,1), send_data(:,2)

  num_of_vgrid = num_of_data

  if (mod(num_of_vgrid,2) == 1) then
    !call jcup_error("moj_interpolation, interpolate_vector", "num_of_vgrid error")
  end if
 
  num_of_vgrid = num_of_vgrid/2
  
  if (present(coef1)) then
    do k = 1, num_of_vgrid
      do i = 1, num_of_operation
        send_grid = send_data_index(i)
        recv_grid = recv_data_index(i)

        !write(6, *) "vector interpolation ", recv_grid, send_grid, coef1(i), coef2(i), send_data(send_grid, 1), send_data(send_grid, 2)

        recv_data(recv_grid, k               ) = recv_data(recv_grid, k               ) &
                                               + send_data(send_grid, k)*coef1(i) &
                                               - send_data(send_grid, k + num_of_vgrid)*coef2(i)
        recv_data(recv_grid, k + num_of_vgrid) = recv_data(recv_grid, k + num_of_vgrid) &
                                               + send_data(send_grid, k)*coef2(i) &
                                               + send_data(send_grid, k + num_of_vgrid)*coef1(i)
      end do
    end do
  else
    do k = 1, num_of_vgrid
      do i = 1, num_of_operation
         send_grid = send_data_index(i)
         recv_grid = recv_data_index(i)

        !write(6, *) "vector interpolation ", recv_grid, send_grid, send_data(send_grid, 1), send_data(send_grid, 2)

         recv_data(recv_grid, k               ) = recv_data(recv_grid, k               ) &
                                                + send_data(send_grid, k)
         recv_data(recv_grid, k + num_of_vgrid) = recv_data(recv_grid, k + num_of_vgrid) &
                                                + send_data(send_grid, k + num_of_vgrid)
      end do
    end do
  end if

  !write(6, *) "moj_interpolate_3d_vector OUT ", size(recv_data, 1), size(recv_data, 2), num_of_data, recv_data(:,1), recv_data(:,2)

end subroutine interpolate_vector

!=======+=========+=========+=========+=========+=========+=========+=========+
!> interpolate data
!! \protected
subroutine interpolation_old(recv_model_name, send_model_name, mapping_tag, send_data, & 
                         recv_data, num_of_data, data_tag)
  implicit none
  character(len=*), intent(IN) :: recv_model_name, send_model_name
  integer, intent(IN) :: mapping_tag
  real(kind=8), intent(IN) :: send_data(:,:)
  real(kind=8), intent(INOUT) :: recv_data(:,:)
  integer, intent(IN) :: num_of_data
  integer, intent(IN) :: data_tag
  logical :: error_flag
  integer :: send_grid, recv_grid, send_coef, recv_coef
  integer :: i, n

  !write(0, *) "moj_interpolation start ", trim(send_model_name)//", "//trim(recv_model_name)
  !write(0, *) mapping_tag, data_tag, size(send_data,1)!, size(send_data,2), num_of_data
  !write(0, *) "SD1 ",send_data(:,1)
  !!write(0, *) "SD2 ",send_data(:,2)

  if (num_of_data <= 0) return ! 20200514 add

  error_flag = .true.

  intpl_ptr => start_ptr

  do while(associated(intpl_ptr))
    if ((trim(intpl_ptr%send_comp) == trim(send_model_name)).and.(trim(intpl_ptr%recv_comp) == trim(recv_model_name))) then
       if (intpl_ptr%mapping_tag == mapping_tag) then
          error_flag = .false.
          exit
       end if
    end if
    intpl_ptr=> intpl_ptr%next_ptr
  end do

  if (error_flag) then
    write(0,*)             "interpolation error!!! send_model="//trim(send_model_name) &
                           //", recv_model="//trim(recv_model_name)
    write(0,'(A,I5,A,I5)') "                        mapping_tag=",mapping_tag,&
                           ", num_of_interpolation=",num_of_interpolation

    write(0,*) associated(start_ptr)

    intpl_ptr => start_ptr
    do while(associated(intpl_ptr))
       write(0,*) "               send_comp="//trim(intpl_ptr%send_comp) &
                                 //", recv_comp="//trim(intpl_ptr%recv_comp)// &
                ", mapping_tag=", intpl_ptr%mapping_tag
       intpl_ptr => intpl_ptr%next_ptr
    end do

    stop
  end if


  if (intpl_ptr%num_of_operation <= 0) return ! 20200514 add

  recv_data(:,:) = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!! exchange_tagが10000より大ならvectorデータ
  if (data_tag > 10000) then ! vector data
      if (associated(intpl_ptr%coef)) then
        call interpolate_vector(send_data, recv_data, num_of_data, &
                                intpl_ptr%num_of_operation, intpl_ptr%send_data_index, intpl_ptr%recv_data_index, &
                                intpl_ptr%coef_vec(:,1), intpl_ptr%coef_vec(:,2))
    
      else
        call interpolate_vector(send_data, recv_data, num_of_data, &
                                intpl_ptr%num_of_operation, intpl_ptr%send_data_index, intpl_ptr%recv_data_index)
      end if

  else ! scalar data

      if (associated(intpl_ptr%coef)) then
        do n = 1, num_of_data
          do i = 1, intpl_ptr%num_of_operation
             send_grid = intpl_ptr%send_data_index(i)
             recv_grid = intpl_ptr%recv_data_index(i)
             send_coef = intpl_ptr%send_coef_index(i)
             recv_coef = intpl_ptr%recv_coef_index(i)
             recv_data(recv_grid, n) = recv_data(recv_grid, n) + send_data(send_grid, n)*intpl_ptr%coef(i) 
             !write(0, *) n, i, recv_grid, send_grid, intpl_ptr%coef(i), recv_data(recv_grid, n), send_data(send_grid, n)
          end do
        end do
      else
        do n = 1, num_of_data
          do i = 1, intpl_ptr%num_of_operation
             send_grid = intpl_ptr%send_data_index(i)
             recv_grid = intpl_ptr%recv_data_index(i)
             send_coef = intpl_ptr%send_coef_index(i)
             recv_coef = intpl_ptr%recv_coef_index(i)
             recv_data(recv_grid, n) = recv_data(recv_grid, n) + send_data(send_grid, n)
             !write(0, *) n, i, recv_data(recv_grid, n), send_data(send_grid, n)
          end do
        end do
      end if

  end if

  !write(0, *) "RD1 ",recv_data(:,1)
  !write(0, *) "RD2 ",recv_data(:,2)

end subroutine interpolation_old

!=======+=========+=========+=========+=========+=========+=========+=========+
!> interpolate data
!! \protected
subroutine interpolation(recv_model_name, send_model_name, mapping_tag, send_data, & 
                         recv_data, num_of_data, data_tag)
  implicit none
  character(len=*), intent(IN) :: recv_model_name, send_model_name
  integer, intent(IN) :: mapping_tag
  real(kind=8), intent(IN) :: send_data(:,:)
  real(kind=8), intent(INOUT) :: recv_data(:,:)
  integer, intent(IN) :: num_of_data
  integer, intent(IN) :: data_tag
  logical :: error_flag
  integer :: send_grid, recv_grid, send_coef, recv_coef
  integer :: i, j, n

  !write(0, *) "moj_interpolation start ", trim(send_model_name)//", "//trim(recv_model_name)
  !write(0, *) mapping_tag, data_tag, size(send_data,1)!, size(send_data,2), num_of_data
  !write(0, *) "SD1 ",send_data(:,1)
  !!write(0, *) "SD2 ",send_data(:,2)

  if (num_of_data <= 0) return ! 20200514 add

  error_flag = .true.

  intpl_ptr => start_ptr

  do while(associated(intpl_ptr))
    if ((trim(intpl_ptr%send_comp) == trim(send_model_name)).and.(trim(intpl_ptr%recv_comp) == trim(recv_model_name))) then
       if (intpl_ptr%mapping_tag == mapping_tag) then
          error_flag = .false.
          exit
       end if
    end if
    intpl_ptr=> intpl_ptr%next_ptr
  end do

  if (error_flag) then
    write(0,*)             "interpolation error!!! send_model="//trim(send_model_name) &
                           //", recv_model="//trim(recv_model_name)
    write(0,'(A,I5,A,I5)') "                        mapping_tag=",mapping_tag,&
                           ", num_of_interpolation=",num_of_interpolation

    write(0,*) associated(start_ptr)

    intpl_ptr => start_ptr
    do while(associated(intpl_ptr))
       write(0,*) "               send_comp="//trim(intpl_ptr%send_comp) &
                                 //", recv_comp="//trim(intpl_ptr%recv_comp)// &
                ", mapping_tag=", intpl_ptr%mapping_tag
       intpl_ptr => intpl_ptr%next_ptr
    end do

    stop
  end if


  if (intpl_ptr%num_of_operation <= 0) return ! 20200514 add

  recv_data(:,:) = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!! exchange_tagが10000より大ならvectorデータ
  if (data_tag > 10000) then ! vector data
      if (associated(intpl_ptr%coef)) then
        call interpolate_vector(send_data, recv_data, num_of_data, &
                                intpl_ptr%num_of_operation, intpl_ptr%send_data_index, intpl_ptr%recv_data_index, &
                                intpl_ptr%coef_vec(:,1), intpl_ptr%coef_vec(:,2))
    
      else
        call interpolate_vector(send_data, recv_data, num_of_data, &
                                intpl_ptr%num_of_operation, intpl_ptr%send_data_index, intpl_ptr%recv_data_index)
      end if

  else ! scalar data

      if (associated(intpl_ptr%coef)) then
        do n = 1, num_of_data
           do j = 1, intpl_ptr%num_of_parallel_intpl
            !$omp parallel do default(none),private(i, send_grid, recv_grid), &
            !$omp shared(j, n, intpl_ptr, send_data, recv_data)
              do i = 1, intpl_ptr%parallel_intpl(j)%num_of_operation
                send_grid = intpl_ptr%parallel_intpl(j)%send_data_index(i)
                recv_grid = intpl_ptr%parallel_intpl(j)%recv_data_index(i)
                recv_data(recv_grid, n) = recv_data(recv_grid, n) + send_data(send_grid, n)*intpl_ptr%parallel_intpl(j)%coef(i) 
                !write(0, *) n, i, recv_grid, send_grid, intpl_ptr%coef(i), recv_data(recv_grid, n), send_data(send_grid, n)
             end do
          end do
        end do
      else
        do n = 1, num_of_data
           do j = 1, intpl_ptr%num_of_parallel_intpl
            !$omp parallel do default(none),private(i, send_grid, recv_grid), &
            !$omp shared(j, n, intpl_ptr, send_data, recv_data)
              do i = 1, intpl_ptr%parallel_intpl(j)%num_of_operation
                 send_grid = intpl_ptr%parallel_intpl(j)%send_data_index(i)
                 recv_grid = intpl_ptr%parallel_intpl(j)%recv_data_index(i)
                 recv_data(recv_grid, n) = recv_data(recv_grid, n) + send_data(send_grid, n)
                 !write(0, *) n, i, recv_data(recv_grid, n), send_data(send_grid, n)
              end do
          end do
        end do
      end if

  end if

  !write(0, *) "RD1 ",recv_data(:,1)
  !write(0, *) "RD2 ",recv_data(:,2)

end subroutine interpolation

!=======+=========+=========+=========+=========+=========+=========+=========+


end module h3ou_intpl

