!====================================================================================================
!> @brief
!> h3open_mp many to one ensemble module
!

module mod_many_to_one
  use mpi
  use mod_ensemble_base
  implicit none
  private

!--------------------------------   public  ----------------------------------!

  public :: set_many_to_one_ensemble  ! subroutine (comp_name, num_of_ensemble)
  
!--------------------------------   private  ---------------------------------!
  
  integer :: ierror
  integer :: errorcode

contains

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_many_to_one_ensemble(comp_name, num_of_ensemble)
  implicit none
  character(len=*), intent(IN) :: comp_name
  integer, intent(IN) :: num_of_ensemble
  integer :: start_rank, next_rank
  integer :: ecomm_id
  

  ecomp%name = trim(comp_name)

  start_rank = 0
  ecomm_id  = 1

  do while(start_rank < global%size)
    call set_ecomp_rank(trim(comp_name), ecomm_id, start_rank, next_rank)
    start_rank = next_rank
    ecomm_id = ecomm_id + 1
  end do

  if (mod(ecomp%size, num_of_ensemble) /= 0) then
     write(0, *) "ERROR!!!, MPI size and ensemble size mismatch"
     call MPI_abort(global%comm, errorcode, ierror)
  end if
  
  call set_local_name(num_of_ensemble)

  call mpi_comm_split(global%comm, ecomp%comm_id, global%rank, ecomp%comm, ierror)
  call mpi_comm_size(ecomp%comm, ecomp%size, ierror)
  call mpi_comm_rank(ecomp%comm, ecomp%rank, ierror)
  call mpi_comm_group(ecomp%comm, ecomp%group, ierror)

  if (couple%comm_id == 1) then ! left side model only
    call mpi_comm_split(global%comm, couple%comm_id, global%rank, couple%comm, ierror)
    call mpi_comm_size(couple%comm, couple%size, ierror)
    call mpi_comm_rank(couple%comm, couple%rank, ierror)
    call mpi_comm_group(couple%comm, couple%group, ierror)
  else
    call set_coupling_flag(.false.)
    call mpi_comm_split(global%comm, couple%comm_id, global%rank, couple%comm, ierror)
  end if
  
  call set_target_name()
  
end subroutine set_many_to_one_ensemble

!=======+=========+=========+=========+=========+=========+=========+=========+

end module mod_many_to_one

