module mod_utils
  use mpi
  private

!--------------------------------   public  ----------------------------------!

  integer, parameter, public :: STR_SHORT = 64
  integer, parameter, public :: STR_MID   = 256
  integer, parameter, public :: STR_LONG  = 1024

  public :: init_utils                   ! subroutine ()
  public :: set_fid                      ! subroutine (fid)
  public :: open_log_file                ! subroutine (basename, my_rank)
  public :: put_log                      ! subroutine (log_str)
  public :: close_log_file               ! subroutine ()
  public :: put_error                    ! subroutine (routine_name, message)
  
!--------------------------------  private  ----------------------------------!

  integer, parameter :: INTVL_SEC = 1
  integer, parameter :: INTVL_MIN = 2
  integer, parameter :: INTVL_HUR = 3
  integer, parameter :: INTVL_DAY = 4
  integer, parameter :: INTVL_MON = 5

  character(len=STR_LONG) :: log_fname ! log file name
  integer :: log_fid = 901             ! log file id
  integer :: log_level = 0
  
  integer :: my_rank_global
  integer :: my_rank_local
  logical :: is_mpi
  integer :: ierror
  integer, parameter :: MIN_FID = 10
  integer, parameter :: MAX_FID = 999

contains

!=======+=========+=========+=========+=========+=========+=========+=========+
!>

subroutine init_utils()
  implicit none

  call mpi_initialized(is_mpi, ierror)

  if (is_mpi) then
    call mpi_comm_rank(MPI_COMM_WORLD, my_rank_global, ierror)
  else
    my_rank_global = -1
  end if

end subroutine init_utils

!=======+=========+=========+=========+=========+=========+=========+=========+
!>

subroutine set_fid(fid)
  implicit none
  integer, intent(INOUT) :: fid
  logical :: op
  
  fid = max(fid, MIN_FID)

  do 
    if (fid > MAX_FID) then
      write(0, *) "[set_fid], fid exceeded MAX_FID"
      stop
    end if
    inquire(unit = fid, OPENED = op)
    if (op) then
       fid = fid + 1
    else
      return
    end if
  end do

end subroutine set_fid

!=======+=========+=========+=========+=========+=========+=========+=========+
!>

subroutine open_log_file(base_name, my_rank, out_log_level)
  implicit none
  character(len=*), intent(IN) :: base_name
  integer, intent(IN) :: my_rank
  integer, intent(IN) :: out_log_level
  
  my_rank_local = my_rank

  log_level = out_log_level
  
  if (log_level <= 0) return
  
  write(log_fname, '(A,A,I0.5)') trim(base_name),".log.PE", my_rank
 
  call set_fid(log_fid)

  open(log_fid, file=trim(log_fname),form = 'formatted', &
           access = 'sequential', action = 'write', err = 200)

  return

200 call put_error('set_log_file','cannot create log file: '//trim(log_fname))

    return
  

end subroutine open_log_file

!=======+=========+=========+=========+=========+=========+=========+=========+
!>

subroutine close_log_file()
  implicit none

  if (log_level <= 0) return
  
  close(log_fid)

end subroutine close_log_file

!=======+=========+=========+=========+=========+=========+=========+=========+
!>

subroutine put_log(log_str, my_log_level)
  implicit none
  character(len=*),  intent(IN) :: log_str
  integer, optional, intent(IN) :: my_log_level
  integer :: out_log_level

  if (present(my_log_level)) then
     out_log_level = my_log_level
  else
     out_log_level = 1
  end if
  
  if (out_log_level <= log_level) write(log_fid, *) trim(log_str)

end subroutine put_log

!=======+=========+=========+=========+=========+=========+=========+=========+
!>

subroutine put_error(routine_name, message)
  implicit none 
  character(len=*),intent(IN) :: routine_name,message

  character(len=STR_LONG) :: message_str
  character(len=STR_LONG) :: istr
  integer :: error_code

  write(istr, *) my_rank_global

  write(message_str, &
      '("!!! error !!! [RANK=",A,"][",A,"] : ",A,", program terminated")') &
      trim(adjustl(istr)),trim(routine_name), trim(message)

  write(0, *) trim(message_str)
  write(log_fid, *) trim(message_str)
  close(log_fid)
  call mpi_abort(MPI_COMM_WORLD, error_code, ierror)
  stop

end subroutine put_error

  

end module mod_utils
