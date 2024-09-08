!====================================================================================================
!> @brief
!> h3ous simple interface for h3-Open-UTIL/MP coupler
!
!Copyright (c) 2024, arakawa@climtech.jp
!All rights reserved.
!
module h3ous
  implicit none
  private
  
!--------------------------------   public  ----------------------------------!

  public :: h3ous_init                      ! subroutine (comp_name)
  public :: h3ous_send                      ! subroutine (target_comp_name, send_data)
  public :: h3ous_recv                      ! subroutine (source_comp_name, recv_data)
  public :: h3ous_end                       ! subroutine ()
  
!--------------------------------   private  ---------------------------------!

  interface h3ous_send
     module procedure h3ous_send_int_scalar
     module procedure h3ous_send_real_scalar
     module procedure h3ous_send_double_scalar
     module procedure h3ous_send_int_array
     module procedure h3ous_send_real_array
     module procedure h3ous_send_double_array
     module procedure h3ous_send_str_array
  end interface h3ous_send

  interface h3ous_recv
     module procedure h3ous_recv_int_scalar
     module procedure h3ous_recv_real_scalar
     module procedure h3ous_recv_double_scalar
     module procedure h3ous_recv_int_array
     module procedure h3ous_recv_real_array
     module procedure h3ous_recv_double_array
     module procedure h3ous_recv_str_array
  end interface h3ous_recv
  
  character(len=64) :: my_name
  integer :: my_rank_global
  integer :: my_size_global
  integer :: ierror
  
contains

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_init(comp_name)
  use mpi
  use h3ou_api, only : h3ou_init
  implicit none
  character(len=*), intent(IN) :: comp_name
  integer :: fid = 128
  integer :: ndata = 1
  integer :: i
  
  write(0, *) "h3ous start, my name =  ", trim(comp_name)  

  call mpi_init(ierror)
  
  call mpi_comm_rank(MPI_COMM_WORLD, my_rank_global, ierror)
  call mpi_comm_size(MPI_COMM_WORLD, my_size_global, ierror)
  
  write(0, *) "h3ous mpi parameters =  ", my_rank_global, my_size_global

  if (my_rank_global == 0) then
        open(unit = fid, file = "coupling.conf", status="replace", action="write")
        write(fid, *) "&h3ou_coupling"
        write(fid, *) '   log_level = "LOUD"'
        write(fid, *) "   debug_mode = .false."
        write(fid, *) "   stop_step = 0"
        write(fid, *) "&end"
        write(fid, *)

        write(fid, *) '&h3ou_var comp_put = "compA" , comp_get = "compB", '
        write(fid, *) '          grid_put = "compA_grid", grid_get = "compB_grid" &end '
        do i = 1, ndata
           write(fid, '(A,I0.0,A,I0.0,A)') ' &h3ou_var var_put = "varA2B', i, '", var_get = "varA2B', i, &
                        '", grid_intpl_tag=1, intvl=1800, lag= 1, layer=1, flag="SNP", &end'
        end do
        write(fid, *) '&h3ou_var comp_put = "compB" , comp_get = "compA", '
        write(fid, *) '          grid_put = "compB_grid", grid_get = "compA_grid" &end '
        do i = 1, ndata
           write(fid, '(A,I0.0,A,I0.0,A)') ' &h3ou_var var_put = "varB2A', i, '", var_get = "varB2A', i, &
                        '", grid_intpl_tag=1, intvl=1800, lag=-1, layer=1, flag="SNP", &end'
        end do
        write(fid, *) '&h3ou_var comp_put = "compC" , comp_get = "compB", '
        write(fid, *) '          grid_put = "compC_grid", grid_get = "compB_grid" &end '
        do i = 1, ndata
           write(fid, '(A,I0.0,A,I0.0,A)') ' &h3ou_var var_put = "varC2B', i, '", var_get = "varC2B', i, &
                        '", grid_intpl_tag=1, intvl=1800, lag= 1, layer=1, flag="SNP", &end'
        end do
        write(fid, *) '&h3ou_var comp_put = "compB" , comp_get = "compC", '
        write(fid, *) '          grid_put = "compB_grid", grid_get = "compC_grid" &end '
        do i = 1, ndata
           write(fid, '(A,I0.0,A,I0.0,A)') ' &h3ou_var var_put = "varB2C', i, '", var_get = "varB2C', i, &
                        '", grid_intpl_tag=1, intvl=1800, lag=-1, layer=1, flag="SNP", &end'
        end do
        close(fid)
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierror)
  
  my_name = trim(comp_name)
  call h3ou_init(my_name, "coupling.conf")
end subroutine h3ous_init
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_end()
  use h3ou_api, only : h3ou_coupling_end
  implicit none

  call mpi_finalize(ierror)

end subroutine h3ous_end

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_send_int_scalar(target_name, val)
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN) :: val
  integer :: send_array(1)

  send_array(1) = val

  call h3ou_send_array(my_name, target_name, send_array)
  
end subroutine h3ous_send_int_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_send_real_scalar(target_name, val)
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=*), intent(IN) :: target_name
  real(kind=4), intent(IN) :: val
  real(kind=4) :: send_array(1)

  send_array(1) = val

  call h3ou_send_array(my_name, target_name, send_array)
  
end subroutine h3ous_send_real_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_send_double_scalar(target_name, val)
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=*), intent(IN) :: target_name
  real(kind=8), intent(IN) :: val
  real(kind=8) :: send_array(1)

  send_array(1) = val

  call h3ou_send_array(my_name, target_name, send_array)
  
end subroutine h3ous_send_double_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_recv_int_scalar(source_name, val)
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=*), intent(IN) :: source_name
  integer, intent(INOUT) :: val
  integer :: recv_array(1)


  call h3ou_recv_array(my_name, source_name, recv_array)
  
  val = recv_array(1)
  
end subroutine h3ous_recv_int_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_recv_real_scalar(source_name, val)
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=*), intent(IN) :: source_name
  real(kind=4), intent(INOUT) :: val
  real(kind=4) :: recv_array(1)


  call h3ou_recv_array(my_name, source_name, recv_array)
  
  val = recv_array(1)
  
end subroutine h3ous_recv_real_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_recv_double_scalar(source_name, val)
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=*), intent(IN) :: source_name
  real(kind=8), intent(INOUT) :: val
  real(kind=8) :: recv_array(1)

  call h3ou_recv_array(my_name, source_name, recv_array)
  
  val = recv_array(1)
  
end subroutine h3ous_recv_double_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_send_int_array(target_name, val)
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN) :: val(:)

  call h3ou_send_array(my_name, target_name, val)
  
end subroutine h3ous_send_int_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_send_real_array(target_name, val)
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=*), intent(IN) :: target_name
  real(kind=4), intent(IN) :: val(:)

  call h3ou_send_array(my_name, target_name, val)
  
end subroutine h3ous_send_real_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_send_double_array(target_name, val)
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=*), intent(IN) :: target_name
  real(kind=8), intent(IN) :: val(:)

  call h3ou_send_array(my_name, target_name, val)
  
end subroutine h3ous_send_double_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_send_str_array(target_name, val)
  use h3ou_api, only : h3ou_send_array
  implicit none
  character(len=*), intent(IN) :: target_name
  character(len=*), intent(IN) :: val

  call h3ou_send_array(my_name, target_name, val)
  
end subroutine h3ous_send_str_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_recv_int_array(source_name, val)
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=*), intent(IN) :: source_name
  integer, intent(INOUT) :: val(:)

  call h3ou_recv_array(my_name, source_name, val)
  
end subroutine h3ous_recv_int_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_recv_real_array(source_name, val)
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=*), intent(IN) :: source_name
  real(kind=4), intent(INOUT) :: val(:)

  call h3ou_recv_array(my_name, source_name, val)
  
end subroutine h3ous_recv_real_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_recv_double_array(source_name, val)
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=*), intent(IN) :: source_name
  real(kind=8), intent(INOUT) :: val(:)

  call h3ou_recv_array(my_name, source_name, val)
  
end subroutine h3ous_recv_double_array

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ous_recv_str_array(source_name, val)
  use h3ou_api, only : h3ou_recv_array
  implicit none
  character(len=*), intent(IN) :: source_name
  character(len=*), intent(INOUT) :: val

  call h3ou_recv_array(my_name, source_name, val)
  
end subroutine h3ous_recv_str_array

!=======+=========+=========+=========+=========+=========+=========+=========+


end module h3ous
