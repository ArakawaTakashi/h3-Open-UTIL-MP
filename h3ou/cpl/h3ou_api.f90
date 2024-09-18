!====================================================================================================
!> @brief
!> h3ou_api jcup interface module 
!
!Copyright (c) 2011, arakawa@rist.jp
!All rights reserved.
!
module h3ou_api
  use mpi
  use jcup_interface, only : jcup_varp_type, jcup_varg_type, &
                             OPERATION_COEF, SEND_COEF, RECV_COEF, &
                             h3ou_inc_calendar => jcup_inc_calendar, &
                             h3ou_get_comp_num_from_name => jcup_get_comp_num_from_name
  use mod_utils, only : STR_SHORT, STR_MID, STR_LONG
  use mod_ensemble_base, only : NAME_LEN, STR_LEN, &
                                NO_ENSEMBLE, ONE_TO_ONE, MANY_TO_ONE, &
                                local, &
                                h3ou_get_my_name => get_my_name, &
                                h3ou_get_target_name => get_target_name, &
                                h3ou_get_couple_id => get_couple_id, &
                                h3ou_is_local_leader => is_local_leader, &
                                h3ou_send_local => send_local, &
                                h3ou_recv_local => recv_local
  use h3ou_intpl, only : h3ou_interpolation => interpolation  
  implicit none
  private

!--------------------------------   public  ----------------------------------!

  public :: jcup_varp_type
  public :: jcup_varg_type
  public :: OPERATION_COEF
  public :: SEND_COEF
  public :: RECV_COEF
  public :: h3ou_init                      ! subroutine (comp_name, config_file_name, num_of_ensemble)
  public :: h3ou_get_my_name               ! character(len=NAME_LEN) function ()
  public :: h3ou_get_target_name           ! character(len=NAME_LEN) function ()
  public :: h3ou_is_coupled                ! logical function()
  public :: h3ou_get_couple_id             ! integer function ()
  public :: h3ou_get_mpi_parameter         ! subroutine (comp_name, comm, group, size, rank)
  public :: h3ou_get_my_comm               ! integer function ()
  public :: h3ou_get_my_rank               ! integer function ()
  public :: h3ou_get_my_size               ! integer function ()
  public :: h3ou_is_local_leader           ! logical function ()
  public :: h3ou_def_grid                  ! subroutine (grid_index, comp_name, grid_name, nz)
  public :: h3ou_end_grid_def              ! subroutine ()
  !public :: h3ou_set_default_configuration ! subroutine (my_comp, send_comp, recv_mode, interval, time_lag, mapping_tag, exchange_tag)
  !public :: h3ou_def_varp                  ! subroutine (data_type_ptr, comp_name, data_name, grid_name, num_of_data)
  !public :: h3ou_def_varg                  ! subroutine (data_type_ptr, 
  !public :: h3ou_end_var_def               ! subroutine ()
  public :: h3ou_get_num_of_put_data       ! integer function ()
  public :: h3ou_get_put_data_name         ! character(len=STR_LEN) function (data_num)
  public :: h3ou_get_num_of_get_data       ! integer function ()
  public :: h3ou_get_get_data_name         ! character(len=STR_LEN) function (data_num)
  public :: h3ou_get_vlayer                ! integer function (data_name)
  public :: h3ou_get_model_id              ! subroutine (comp_name, comp_id)
  public :: h3ou_get_comp_num_from_name    ! integer function (comp_name)
  public :: h3ou_is_my_component           ! logical funcion (comp_id)
  public :: h3ou_set_interpolation_table   ! subroutine (my_comp, send_comp, send_grid, recv_comp, recv_grid, tag,
                                           !             send_index, recv_indec, coef, coef_cos, coef_sin)
  !public :: h3ou_set_mapping_table         ! subroutine (my_comp, send_comp, send_grid, recv_comp, recv_grid, tag, send_index, recv_index)  
  !public :: h3ou_set_local_coef            ! subroutine (my_comp, send_comp, mapping_tag, global_coef, local_coef, coef_type)
  !public :: h3ou_send_coef                 ! subroutine (my_comp, recv_comp, coef)
  !public :: h3ou_recv_coef                 !    subroutine (my_comp, send_comp, coef)
  !public :: h3ou_get_local_operation_index ! subroutine (
  public :: h3ou_get_num_of_send_grid      ! subroutine (recv_comp, send_comp, grid_tag, num_of_grid)
  public :: h3ou_get_num_of_recv_grid      ! subroutine (recv_comp, send_comp, grid_tag, num_of_grid)
  public :: h3ou_init_time                 ! subroutine (time_array)
  public :: h3ou_set_time                  ! subroutine (comp_name, time_array, delta_t, is_exchange)
  public :: h3ou_inc_calendar              ! subroutine (time_array, delta_t)
  public :: h3ou_inc_time                  ! subroutine (comp_name, time_array)
  public :: h3ou_send_data_immediately     ! subroutine (send_comp, recv_comp, time_lag)
  public :: h3ou_recv_data_immediately     ! subroutine (send_comp, recv_comp)
  public :: h3ou_put_data                  ! subroutine (varp_type, data, data_vector)
  public :: h3ou_put_data_1d               ! subroutine (varp_type, data, data_vector)
  public :: h3ou_put_data_25d              ! subroutine (varp_type, data, data_vector)
  public :: h3ou_get_data                  ! subroutine (varg_type, data, data_vector, is_recv_ok)
  public :: h3ou_get_data_1d               ! subroutine (varg_type, data, data_vector, is_recv_ok)
  public :: h3ou_get_data_25d              ! subroutine (varg_type, data, data_vector, is_recv_ok)
  public :: h3ou_send_local                ! subroutine
  public :: h3ou_recv_local                ! subroutine
  public :: h3ou_send_array                ! subroutine (my_comp_name, recv_comp_name, array)
  public :: h3ou_recv_array                ! subroutine (my_comp_name, send_comp_name, array)
  public :: h3ou_coupling_end              ! subroutine ()
  public :: h3ou_interpolation

  !---------------------------------------------------------------------------
  public :: h3ou_end                       ! subroutine ()
  public :: h3ou_bcast_global_string       ! subroutine (source_pe, str_data)
  public :: h3ou_bcast_global_int          ! subroutine (source_pe, int_data)
  public :: h3ou_bcast_global_real         ! subroutine (source_pe, real_data)
  public :: h3ou_bcast_global_double       ! subroutine (source_pe, double_data)
  public :: h3ou_bcast_global              ! subroutine (source_pe, data)
  public :: h3ou_bcast_model_int           ! subroutine (source_name, target_name, data)
  public :: h3ou_bcast_model_real          ! subroutine (source_name, target_name, data)
  public :: h3ou_bcast_model_double        ! subroutine (source_name, target_name, data)
  public :: h3ou_bcast_model               ! subroutine (source_name, target_name, data)
  public :: h3ou_send_model_int            ! subroutine (target_name, target_pe, data)
  public :: h3ou_send_model_real           ! subroutine (target_name, target_pe, data)
  public :: h3ou_send_model_double         ! subroutine (target_name, target_pe, data)
  public :: h3ou_recv_model_int            ! subroutine (source_name, source_pe, data)
  public :: h3ou_recv_model_real           ! subroutine (source_name, source_pe, data)
  public :: h3ou_recv_model_double         ! subroutine (source_name, source_pe, data)
  public :: h3ou_send_model                ! subroutine (target_name, target_pe, data)
  public :: h3ou_recv_model                ! subroutine (source_name, source_pe, data)
  public :: h3ou_bcast_local_int           ! subroutine (source_pe, data)
  public :: h3ou_bcast_local_real          ! subroutine (source_pe, data)
  public :: h3ou_bcast_local_double        ! subroutine (source_pe, data)
  public :: h3ou_bcast_local               ! subroutine (source_pe, data)
  public :: h3ou_send_local_int            ! subroutine (target_pe, data)
  public :: h3ou_send_local_real           ! subroutine (target_pe, data)
  public :: h3ou_send_local_double         ! subroutine (target_pe, data)
  public :: h3ou_recv_local_int            ! subroutine (source_pe, data)
  public :: h3ou_recv_local_real           ! subroutine (source_pe, data)
  public :: h3ou_recv_local_double         ! subroutine (source_pe, data)
  public :: h3ou_send_local_all            ! subroutine (target_pe, data)
  public :: h3ou_recv_local_all            ! subroutine (target_pe, data)
  public :: h3ou_send                      ! subroutine (target_comp, val)
  public :: h3ou_recv                      ! subroutine (source_comp, val)
  public :: h3ou_isend_model_int           ! subroutine (target_name, target_pe, data)
  public :: h3ou_isend_model_real          ! subroutine (target_name, target_pe, data)
  public :: h3ou_isend_model_double        ! subroutine (target_name, target_pe, data)
  public :: h3ou_irecv_model_int           ! subroutine (source_name, source_pe, data)
  public :: h3ou_irecv_model_real          ! subroutine (source_name, source_pe, data)
  public :: h3ou_irecv_model_double        ! subroutine (source_name, source_pe, data)
  public :: h3ou_isend_waitall             ! subroutine ()
  public :: h3ou_irecv_waitall             ! subroutine ()
  
!--------------------------------   private  ---------------------------------!
  
  interface h3ou_put_data
     module procedure h3ou_put_data_1d, h3ou_put_data_25d
  end interface h3ou_put_data

  interface h3ou_get_data
     module procedure h3ou_get_data_1d, h3ou_get_data_25d
  end interface h3ou_get_data

  interface h3ou_bcast_global
     module procedure h3ou_bcast_global_string
     module procedure h3ou_bcast_global_int
     module procedure h3ou_bcast_global_real
     module procedure h3ou_bcast_global_double
  end interface h3ou_bcast_global

  interface h3ou_bcast_model
     module procedure h3ou_bcast_model_int
     module procedure h3ou_bcast_model_real
     module procedure h3ou_bcast_model_double
  end interface h3ou_bcast_model
  
  interface h3ou_send_array
     module procedure h3ou_send_array_str
     module procedure h3ou_send_array_int
     module procedure h3ou_send_array_real
     module procedure h3ou_send_array_dbl
     module procedure h3ou_send_array_str_2
     module procedure h3ou_send_array_int_2
     module procedure h3ou_send_array_real_2
     module procedure h3ou_send_array_dbl_2
  end interface h3ou_send_array

  interface h3ou_recv_array
     module procedure h3ou_recv_array_str
     module procedure h3ou_recv_array_int
     module procedure h3ou_recv_array_real
     module procedure h3ou_recv_array_dbl
     module procedure h3ou_recv_array_str_2
     module procedure h3ou_recv_array_int_2
     module procedure h3ou_recv_array_real_2
     module procedure h3ou_recv_array_dbl_2
  end interface h3ou_recv_array

  interface h3ou_send
     module procedure h3ou_send_int_scalar
     module procedure h3ou_send_real_scalar
     module procedure h3ou_send_double_scalar
     module procedure h3ou_send_array_str
     module procedure h3ou_send_array_int
     module procedure h3ou_send_array_real
     module procedure h3ou_send_array_dbl
  end interface h3ou_send

  interface h3ou_recv
     module procedure h3ou_recv_int_scalar
     module procedure h3ou_recv_real_scalar
     module procedure h3ou_recv_double_scalar
     module procedure h3ou_recv_array_str
     module procedure h3ou_recv_array_int
     module procedure h3ou_recv_array_real
     module procedure h3ou_recv_array_dbl
  end interface h3ou_recv

  interface h3ou_send_model
     module procedure h3ou_send_model_int
     module procedure h3ou_send_model_real
     module procedure h3ou_send_model_double
  end interface h3ou_send_model
  
  interface h3ou_recv_model
     module procedure h3ou_recv_model_int
     module procedure h3ou_recv_model_real
     module procedure h3ou_recv_model_double
  end interface h3ou_recv_model

  interface h3ou_bcast_local
     module procedure h3ou_bcast_local_int
     module procedure h3ou_bcast_local_real
     module procedure h3ou_bcast_local_double
  end interface h3ou_bcast_local

  interface h3ou_send_local_all
     module procedure h3ou_send_local_int
     module procedure h3ou_send_local_real
     module procedure h3ou_send_local_double
  end interface h3ou_send_local_all
  
  interface h3ou_recv_local_all
     module procedure h3ou_recv_local_int
     module procedure h3ou_recv_local_real
     module procedure h3ou_recv_local_double
  end interface h3ou_recv_local_all
  
     
  character(len=STR_LONG) :: log_str
    
  integer :: ierror
  integer :: errorcode
  integer :: stop_step = 0
  
  logical :: is_coupled

  character(len=STR_SHORT) :: my_name ! my component name 
  integer                  :: my_id   ! id of my component

  integer :: local_comm, local_group, local_size, local_rank
  
  logical :: is_full_coupling = .true. ! full coupling mode or simple exchange mode

contains

!=======+=========+=========+=========+=========+=========+=========+=========+
!> @breaf
!> set component name
!> @param[in] component_name name of component
subroutine h3ou_init(comp_name, config_file_name, nensemble)
  use mod_utils, only : open_log_file, put_log, close_log_file
  use mod_namelist, only : read_coupler_config, get_log_level, is_debug_mode, read_namelist, get_num_of_configuration, get_stop_step
  use mod_comp, only : init_comp
  use mod_ensemble_base, only : NO_ENSEMBLE, &
                                eb_init, eb_end_init, get_ensemble_type, get_coupling_flag, get_couple_id, &
                                get_my_name, get_target_name, get_local_mpi_rank, get_local_mpi_param
  use mod_one_to_one, only : set_one_to_one_ensemble
  use mod_many_to_one, only : set_many_to_one_ensemble
  use h3ou_intpl, only : init_interpolation
  use jcup_interface, only : jcup_set_new_comp, jcup_initialize, jcup_log, jcup_get_model_id
  use palmtime, only : palm_TimeInit
  implicit none
  character(len=*), intent(IN) :: comp_name
  character(len=*), intent(IN) :: config_file_name
  integer, optional, intent(IN) :: nensemble
  character(len=NAME_LEN) :: e_comp_name ! ensemble component name
  integer :: num_of_ensemble
  logical :: init_flag
  integer :: start_rank, next_rank
  integer :: glocal_id
  integer :: i
  character(len=STR_LONG) :: log_str

  if (present(nensemble)) then
     num_of_ensemble = nensemble
  else
     num_of_ensemble = 1
  end if

  call eb_init(comp_name, num_of_ensemble)
  
  is_coupled = .true.

  select case(get_ensemble_type())
    case (NO_ENSEMBLE, ONE_TO_ONE)
      call set_one_to_one_ensemble(comp_name, num_of_ensemble)
    case (MANY_TO_ONE)
      call set_many_to_one_ensemble(comp_name, num_of_ensemble)
  end select

  is_coupled = get_coupling_flag()

  call read_coupler_config(trim(config_file_name))

  
  call eb_end_init(is_debug_mode())

  call get_local_mpi_param(local_comm, local_group, local_size, local_rank)

  call palm_TimeInit(trim(comp_name), local_comm)
  
  my_name = trim(get_my_name())
  
  call open_log_file("h3ou."//trim(get_my_name()), get_local_mpi_rank(), get_log_level())

  stop_step = get_stop_step()

  if (stop_step == 1) then
    call put_log("[stop_step = 1], h3ou_init:IN   OK")
    call close_log_file()
    call mpi_finalize(ierror)
    write(0, *) "11111 --- Coupler stopped at the beginning of h3ou_init --- 11111"
    stop
  end if

  call put_log("=============================================================================")
  call put_log("===============       h3-Open-UTIL/MP coupling log       ====================")
  call put_log("=============================================================================")
  call put_log(">>>>>>>>>>>>>>> h3ou_init IN")

  if (get_ensemble_type() /= NO_ENSEMBLE) then

    call put_log("--------------- ensemble information")
    write(log_str, '("Ensemble coupling : my name = ", A10, ", target name = ", A10)') trim(get_my_name()), trim(get_target_name())

    call put_log(trim(log_str))
    call jcup_log("h3ou_init", trim(log_str), 2)

    write(log_str, '("ensemble type = ", I1, ", coupling flag = ", L)') get_ensemble_type(), get_coupling_flag()
 
    call put_log(trim(log_str))
    call jcup_log("h3ou_init", trim(log_str), 2)
  end if
 
  call read_namelist(comp_name)

  if (get_num_of_configuration() == 0) then
     is_full_coupling = .false.
     call put_log(" coupling mode : simple exchange mode")
  else
     is_full_coupling = .true.
     call put_log(" coupling mode : full coupling mode")
  end if
  
  call init_comp(get_my_name())

  call jcup_get_model_id(trim(my_name), my_id)
  
  if (is_full_coupling) call init_interpolation()

  call put_log("<<<<<<<<<<<<<<< h3ou_init OUT")

  if (stop_step == 2) then
    call put_log("[stop_step = 2], h3ou_init:OUT  OK. Breakpoint check run completed.")
    call close_log_file() 
    call mpi_finalize(ierror)
    write(0, *) "22222 --- Coupler stopped at the end of h3ou_init --- 22222"
    stop 
  end if

  return

end subroutine h3ou_init

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_coupling_end(itime, is_call_finalize)
  use jcup_interface, only : jcup_coupling_end, jcup_log
  use palmtime, only : palm_TimeFinalize
  use mod_utils, only : put_log, close_log_file
  implicit none
  integer, optional, intent(IN) :: itime(:)
  logical, optional, intent(IN) :: is_call_finalize
  integer :: ierr

  call put_log(">>>>>>>>>>>>>>> h3ou_coupling_end IN")
  call palm_TimeFinalize()  ! not output TM* files

  if (is_full_coupling) then
    if (is_coupled) then
      call jcup_coupling_end(itime, is_call_finalize)
    else
       if (present(is_call_finalize)) then
          if (is_call_finalize) then
             call jcup_log("h3ou_coupling_end", "ensemble coupling finish")
             call mpi_finalize(ierr)
          end if
      end if
    end if
  else
     call mpi_finalize(ierr)
  end if
  call put_log("<<<<<<<<<<<<<<< h3ou_coupling_end OUT")

  call close_log_file()

end subroutine h3ou_coupling_end

!=======+=========+=========+=========+=========+=========+=========+=========+

function h3ou_is_coupled() result(res)
  implicit none
  logical :: res

  res = is_coupled

end function h3ou_is_coupled

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_get_mpi_parameter(comp_name, comm, group, size, rank)
  use mod_ensemble_base, only : get_local_mpi_param
  implicit none
  character(len=*), intent(IN) :: comp_name
  integer, intent(OUT) :: comm, group, size, rank

  call get_local_mpi_param(comm, group, size, rank)

end subroutine h3ou_get_mpi_parameter

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function h3ou_get_my_comm()
  implicit none

  h3ou_get_my_comm = local_comm

end function h3ou_get_my_comm

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function h3ou_get_my_rank()
  implicit none

  h3ou_get_my_rank = local_rank

end function h3ou_get_my_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function h3ou_get_my_size()
  implicit none

  h3ou_get_my_size = local_size

end function h3ou_get_my_size

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_def_grid(grid_index, comp_name, grid_name, nz)
  use mod_utils, only : put_log, close_log_file
  use jcup_interface, only : jcup_def_grid
  use palmtime, only : palm_TimeStart, palm_TimeEnd
  use mod_ensemble_base, only : get_my_name
  use mod_utils, only : put_log
  implicit none
  integer, intent(IN) :: grid_index(:)
  character(len=*), intent(IN) :: comp_name
  character(len=*), intent(IN) :: grid_name
  integer, intent(IN)  :: nz
  character(len=STR_LONG) :: log_str
  
  call put_log("")
  call put_log(">>>>>>>>>>>>>>> h3ou_def_grid IN")
  call palm_TimeStart("def_grid")

 if (stop_step == 3) then
    call put_log("[stop_step = 3], h3ou_def_grid:IN   OK")
    call close_log_file()
    call mpi_finalize(ierror)
    write(0, *) "33333 --- Coupler stopped at the beginning of h3ou_def_grid --- 33333"
    stop 
  end if

  if (is_coupled) then
     call jcup_def_grid(grid_index, trim(get_my_name()), grid_name, nz)
  end if

  call palm_TimeEnd("def_grid")

  call put_log("--------------- grid information")
  write(log_str, '(A,I0,A,I0,A,I0,A,I0)') " grid_name = "//trim(grid_name)//", grid_array_size = ", size(grid_index), &
                                     ", grid_min = ",minval(grid_index), ", grid_max = ", maxval(grid_index),", nz = ", nz
  call put_log(trim(log_str))
  call put_log("<<<<<<<<<<<<<<< h3ou_def_grid OUT")

  if (stop_step == 4) then
    call put_log("[stop_step = 4], h3ou_def_grid:OUT  OK. Breakpoint check run completed.")
    call close_log_file()
    call mpi_finalize(ierror)
    write(0, *) "44444 --- Coupler stopped at the end of h3ou_def_grid --- 44444"
    stop 
  end if

end subroutine h3ou_def_grid

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_end_grid_def()
  use mod_utils, only : put_log, close_log_file
  use jcup_interface, only : jcup_end_grid_def
  use mod_comp, only : comp_type, varp_type, varg_type, get_num_of_comp, get_comp_ptr, &
                       get_varp_ptr, get_num_of_varp, get_varg_ptr, get_num_of_varg
  use mod_utils, only : put_log
  use palmtime, only : palm_TimeStart, palm_TimeEnd
  implicit none
  type(comp_type), pointer :: cptr
  type(varp_type), pointer :: varp_ptr
  type(varg_type), pointer :: varg_ptr
  character(len=STR_LONG) :: log_str
  integer :: i, j

  if (.not.is_coupled) return

  call put_log("")
  call put_log(">>>>>>>>>>>>>>> h3ou_end_grid_def IN")

  if (stop_step == 5) then
    call put_log("[stop_step = 5], h3ou_end_grid_def:IN   OK")
    call close_log_file()
    call mpi_finalize(ierror)
    write(0, *) "55555 --- Coupler stopped at the beginning of h3ou_end_grid_def --- 55555"
    stop 
  end if


  call palm_TimeStart("end_grid_def")

  call jcup_end_grid_def()

  call put_log("--------------- exchange data information")
  do i = 1, get_num_of_comp()
     cptr => get_comp_ptr(i)
     do j = 1, get_num_of_varp(cptr)
        varp_ptr => get_varp_ptr(cptr%pdata_ptr, j)

        call h3ou_def_varp(varp_ptr%jcup_varp_ptr, trim(cptr%name), trim(varp_ptr%name), &
                          trim(varp_ptr%grid_name), varp_ptr%num_of_layer)
        write(log_str, '(A, I0)') " send_data_name = "//trim(varp_ptr%name)//", num_of_layer = ", varp_ptr%num_of_layer
        call put_log(trim(log_str))
     end do
     do j = 1, get_num_of_varg(cptr)
        varg_ptr => get_varg_ptr(cptr%gdata_ptr, j)
        call h3ou_def_varg(varg_ptr%jcup_varg_ptr, trim(cptr%name), trim(varg_ptr%name), &
                          trim(varg_ptr%my_grid_name), varg_ptr%num_of_layer, &
                          trim(varg_ptr%send_comp), trim(varg_ptr%send_data), &
                          varg_ptr%flag, varg_ptr%intvl, varg_ptr%time_lag, 1, varg_ptr%data_tag)
        write(log_str, '(A,I0,A,I0,A,I0)') " recv_data_name = "//trim(varg_ptr%name)//", num_of_layer = ", varg_ptr%num_of_layer, &
                                           ", exchange_interval = ", varg_ptr%intvl, ", send_data_name = "//trim(varg_ptr%send_data)
        call put_log(trim(log_str))
     !send_model_name, send_data_name, recv_mode, interval, time_lag, &
     !mapping_tag, exchange_tag)
      end do
     
  end do

  call h3ou_end_var_def()
  
  call palm_TimeEnd("end_grid_def")
  call put_log("<<<<<<<<<<<<<<< h3ou_end_grid_def OUT")

  if (stop_step == 6) then
    call put_log("[stop_step = 6], h3ou_end_grid_def:OUT  OK. Breakpoint check run completed.")
    call close_log_file()
    call mpi_finalize(ierror)
    write(0, *) "66666 --- Coupler stopped at the end of h3ou_end_grid_def --- 66666"
    stop 
  end if

end subroutine h3ou_end_grid_def


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_set_default_configuration(my_comp, send_comp, recv_mode, interval, time_lag, mapping_tag, exchange_tag)
  use jcup_interface, only : jcup_set_default_configuration
  use mod_ensemble_base, only : get_my_name
  implicit none
  character(len=*), intent(IN) :: my_comp, send_comp
  character(len=3), optional, intent(IN) :: recv_mode
  integer, optional, intent(IN) :: interval
  integer, optional, intent(IN) :: time_lag
  integer, optional, intent(IN) :: mapping_tag
  integer, optional, intent(IN) :: exchange_tag

  if (is_coupled) &
    call jcup_set_default_configuration(trim(get_my_name()), send_comp, recv_mode, interval, time_lag, mapping_tag, exchange_tag)

end subroutine h3ou_set_default_configuration

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_def_varp(data_type_ptr, comp_name, data_name, grid_name, num_of_data)
  use jcup_interface, only : jcup_def_varp
  implicit none
  type(jcup_varp_type), pointer :: data_type_ptr
  character(len=*), intent(IN) :: comp_name
  character(len=*), intent(IN) :: data_name
  character(len=*), intent(IN) :: grid_name
  integer, optional, intent(IN) :: num_of_data

  if (is_coupled) call jcup_def_varp(data_type_ptr, comp_name, data_name, grid_name, num_of_data)
  
end subroutine h3ou_def_varp

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_def_varg(data_type_ptr, comp_name, data_name, grid_name, num_of_data, &
      send_model_name, send_data_name, recv_mode, interval, time_lag, &
     mapping_tag, exchange_tag)
  use jcup_interface, only : jcup_def_varg
  implicit none
  type(jcup_varg_type), pointer :: data_type_ptr
  character(len=*), intent(IN) :: comp_name
  character(len=*), intent(IN) :: data_name
  character(len=*), intent(IN) :: grid_name
  integer, optional, intent(IN) :: num_of_data
  character(len=*), intent(IN) :: send_model_name
  character(len=*), intent(IN) :: send_data_name
  character(len=3), optional, intent(IN) :: recv_mode
  integer, optional, intent(IN) :: interval
  integer, optional, intent(IN) :: time_lag
  integer, optional, intent(IN) :: mapping_tag
  integer, optional, intent(IN) :: exchange_tag

  if (is_coupled) call jcup_def_varg(data_type_ptr, comp_name, data_name, grid_name, num_of_data, &
       send_model_name, send_data_name, recv_mode, interval, time_lag, &
       mapping_tag, exchange_tag)

end subroutine h3ou_def_varg

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_end_var_def()
  use jcup_interface, only : jcup_end_var_def
  implicit none

  if (is_coupled) call jcup_end_var_def()
  
end subroutine h3ou_end_var_def

!=======+=========+=========+=========+=========+=========+=========+=========+

function h3ou_get_num_of_put_data() result(res)
  use mod_comp, only : comp_type, get_comp_ptr
  implicit none
  type(comp_type), pointer :: comp_ptr
  integer :: res
  
  comp_ptr => get_comp_ptr(trim(my_name))

  res = comp_ptr%num_of_put_data

end function h3ou_get_num_of_put_data

!=======+=========+=========+=========+=========+=========+=========+=========+

function h3ou_get_put_data_name(data_num) result(res)
  use mod_comp, only : comp_type, varp_type, get_comp_ptr, get_varp_ptr 
  implicit none
  integer, intent(IN) :: data_num
  type(comp_type), pointer :: comp_ptr
  type(varp_type), pointer :: varp_ptr
  character(len=STR_SHORT) :: res
  
  comp_ptr => get_comp_ptr(trim(my_name))
  varp_ptr => get_varp_ptr(comp_ptr%pdata_ptr, data_num)
  res = trim(varp_ptr%name)

end function h3ou_get_put_data_name

!=======+=========+=========+=========+=========+=========+=========+=========+

function h3ou_get_num_of_get_data() result(res)
  use mod_comp, only : comp_type, get_comp_ptr
  implicit none
  type(comp_type), pointer :: comp_ptr
  integer :: res
  
  comp_ptr => get_comp_ptr(trim(my_name))

  res = comp_ptr%num_of_get_data

end function h3ou_get_num_of_get_data

!=======+=========+=========+=========+=========+=========+=========+=========+

function h3ou_get_get_data_name(data_num) result(res)
  use mod_comp, only : comp_type, varg_type, get_comp_ptr, get_varg_ptr 
  implicit none
  integer, intent(IN) :: data_num
  type(comp_type), pointer :: comp_ptr
  type(varg_type), pointer :: varg_ptr
  character(len=STR_SHORT) :: res
  
  comp_ptr => get_comp_ptr(trim(my_name))
  varg_ptr => get_varg_ptr(comp_ptr%gdata_ptr, data_num)
  res = trim(varg_ptr%name)

end function h3ou_get_get_data_name

!=======+=========+=========+=========+=========+=========+=========+=========+

function h3ou_get_vlayer(data_name) result(res)
  use mod_comp, only : comp_type, varp_type, varg_type, get_comp_ptr, get_varp_ptr, get_varg_ptr 
  use mod_utils, only : put_error
  implicit none
  character(len=*), intent(IN) :: data_name
  type(comp_type), pointer :: comp_ptr
  type(varp_type), pointer :: varp_ptr
  type(varg_type), pointer :: varg_ptr
  integer :: res
  
  comp_ptr => get_comp_ptr(trim(my_name))

  varp_ptr => get_varp_ptr(comp_ptr%pdata_ptr, trim(data_name))
  if (associated(varp_ptr)) then
     res = varp_ptr%num_of_layer
     return
  end if
  
  varg_ptr => get_varg_ptr(comp_ptr%gdata_ptr, trim(data_name))
  if (associated(varg_ptr)) then
     res = varg_ptr%num_of_layer
     return
  end if

  call put_error("h3ou_get_vlayer", "Cannot find the data : "//trim(data_name))

end function h3ou_get_vlayer


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_get_model_id(comp_name, comp_id)
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: comp_name
  integer, intent(OUT) :: comp_id
  
  if (is_coupled) then
     call jcup_get_model_id(comp_name, comp_id)
  else
     comp_id = 0
  end if
  
end subroutine h3ou_get_model_id


!=======+=========+=========+=========+=========+=========+=========+=========+

function h3ou_is_my_component(comp_id) result(res)
  use jcup_interface, only : jcup_is_my_component
  implicit none
  integer, intent(IN) :: comp_id
  logical :: res

  if (is_coupled) then
     res = jcup_is_my_component(comp_id)
  else
     res = .false.
  end if
  
end function h3ou_is_my_component

!=======+=========+=========+=========+=========+=========+=========+=========+
!>
subroutine h3ou_set_interpolation_table(my_name, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name, &
                                       mapping_tag, send_index, recv_index, coef, coef_cos, coef_sin)
  use jcup_interface, only : jcup_set_mapping_table, jcup_send_array, jcup_recv_array, jcup_error
  use h3ou_intpl, only : set_interpolation
  use palmtime, only : palm_TimeStart, palm_TimeEnd
  use mod_ensemble_base, only : NO_ENSEMBLE, get_ensemble_type, get_my_name, get_target_name
  use mod_utils, only : put_log
  implicit none
  character(len=*), intent(IN) :: my_name
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer, intent(IN) :: mapping_tag
  integer, optional, intent(IN) :: send_index(:), recv_index(:)
  real(kind=8), optional, intent(IN) :: coef(:)
  real(kind=8), optional, intent(IN) :: coef_cos(:), coef_sin(:)
  real(kind=8), allocatable :: my_coef(:), my_coef_cos(:), my_coef_sin(:)
  integer :: int_array(2)
  character(len=NAME_LEN) :: my_ensemble_name, send_ensemble_name, recv_ensemble_name
  character(len=STR_LONG) :: log_str
  
  if (.not.is_coupled) return
  
  call put_log("")
  call put_log(">>>>>>>>>>>>>>> h3ou_set_interpolation_table IN")
  call put_log("--------------- interpolation table information")

  call palm_TimeStart("set_interpolation_table")

  if (get_ensemble_type() /= NO_ENSEMBLE) then
     if (trim(my_name) == trim(send_comp_name)) then
       my_ensemble_name   = trim(get_my_name())
       send_ensemble_name = trim(get_my_name())
       recv_ensemble_name = trim(get_target_name())
     else
       my_ensemble_name   = trim(get_my_name())
       send_ensemble_name = trim(get_target_name())
       recv_ensemble_name = trim(get_my_name())
     end if
  else
     my_ensemble_name   = trim(my_name)
     send_ensemble_name = trim(send_comp_name)
     recv_ensemble_name = trim(recv_comp_name)
  end if

  write(log_str, '(A)') " send_comp_name = "//trim(send_ensemble_name)//", send_grid_name = "//trim(send_grid_name)
  call put_log(trim(log_str))
  write(log_str, '(A)') " recv_comp_name = "//trim(recv_ensemble_name)//", recv_grid_name = "//trim(recv_grid_name)
  call put_log(trim(log_str))

  if (present(send_index)) then
    write(log_str, '(A,I0,A,I0,A,I0,A,I0,A,I0)') " index_size = ", size(send_index), &
                                                 ", send_min = ", minval(send_index),", send_max = ", maxval(send_index), & 
                                                 ", recv_min = ", minval(recv_index),", recv_max = ", maxval(recv_index)
    if (h3ou_is_local_leader()) call put_log(trim(log_str))
  end if
 
  
  if (present(send_index)) then
    call jcup_set_mapping_table(my_ensemble_name, send_ensemble_name, send_grid_name, recv_ensemble_name, recv_grid_name, &
                                mapping_tag, send_index, recv_index)
  else
    call jcup_set_mapping_table(my_ensemble_name, send_ensemble_name, send_grid_name, recv_ensemble_name, recv_grid_name, &
                                mapping_tag)
  end if

  
  int_array(:) = 0
  if (trim(my_name) == trim(send_comp_name)) then ! my comp == send comp

    if (present(coef)) then
       int_array(1) = 1 ; int_array(2) = size(coef)
    end if
    if (present(coef_cos)) then
       int_array(1) = 3
       int_array(2) = size(coef_cos)
    end if

    call jcup_send_array(my_ensemble_name, recv_ensemble_name, int_array) ! send the number of coef and the arry size
 
    if (present(coef)) then
      call jcup_send_array(my_ensemble_name, recv_ensemble_name, coef)
    end if

    if (present(coef_cos)) then
       call jcup_send_array(my_ensemble_name, recv_ensemble_name, coef_cos)
       call jcup_send_array(my_ensemble_name, recv_ensemble_name, coef_sin)
    end if

  else ! my comp == recv comp

    call jcup_recv_array(my_ensemble_name, send_ensemble_name, int_array)

    if (int_array(1) >= 1) then ! exist coef

       if (present(coef).or.present(coef_cos)) then
         call jcup_error("moj_set_interpolation_table", "coef or coef_cos double definition")
       end if

       allocate(my_coef(int_array(2)))
       call jcup_recv_array(my_ensemble_name, send_ensemble_name, my_coef, bcast_flag = .false.)
   end if

    if (int_array(1) == 3) then ! coef_cos + coef_sin     
       if (present(coef).or.present(coef_cos)) then
         call jcup_error("moj_set_interpolation_table", "coef or coef_cos double definition")
       end if
       allocate(my_coef_cos(int_array(2)), my_coef_sin(int_array(2)))
       call jcup_recv_array(my_ensemble_name, send_ensemble_name, my_coef_cos, bcast_flag = .false.)
       call jcup_recv_array(my_ensemble_name, send_ensemble_name, my_coef_sin, bcast_flag = .false.)
    end if

    if (present(coef)) then
       int_array(1) = 1 ; int_array(2) = size(coef)
       allocate(my_coef(int_array(2)))
       my_coef(:) = coef(:)
    end if
 
   if (present(coef_cos)) then
       int_array(1) = 3 ; int_array(2) = size(coef_cos)
       allocate(my_coef_cos(int_array(2)), my_coef_sin(int_array(2)))
       my_coef_cos(:) = coef_cos(:)
       my_coef_sin(:) = coef_sin(:)
    end if

    select case(int_array(1))
    case(0) ! no coef
      call set_interpolation(send_ensemble_name, send_grid_name, recv_ensemble_name, recv_grid_name, mapping_tag)
    case(1) ! coef only 
      call set_interpolation(send_ensemble_name, send_grid_name, recv_ensemble_name, recv_grid_name, mapping_tag, my_coef)
      deallocate(my_coef)
    case(3) ! coef + coef_cos + coef_sin
      call set_interpolation(send_ensemble_name, send_grid_name, recv_ensemble_name, recv_grid_name, mapping_tag, &
                             my_coef, my_coef_cos, my_coef_sin)
      deallocate(my_coef, my_coef_cos, my_coef_sin)

    end select

  end if
   
  call palm_TimeEnd("set_interpolation_table")

  call put_log("<<<<<<<<<<<<<<< h3ou_set_interpolatoin_table OUT")

end subroutine h3ou_set_interpolation_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_set_mapping_table(my_model_name, &
                                 send_model_name, send_grid_name, &
                                 recv_model_name,  recv_grid_name, mapping_tag, &
                                 send_grid, recv_grid)
  use jcup_interface, only : jcup_set_mapping_table
  implicit none
  character(len=*), intent(IN)  :: my_model_name
  character(len=*), intent(IN)  :: send_model_name, send_grid_name
  character(len=*), intent(IN)  :: recv_model_name, recv_grid_name
  integer, intent(IN)           :: mapping_tag
  integer, intent(IN), optional :: send_grid(:), recv_grid(:)

  if (is_coupled) then
    call jcup_set_mapping_table(my_model_name, send_model_name, send_grid_name, &
                                               recv_model_name, recv_grid_name, &
                                               mapping_tag, send_grid, recv_grid)
 end if

end subroutine h3ou_set_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========

subroutine h3ou_set_local_coef(my_comp, send_comp, mapping_tag, global_coef, &
                              local_coef, coef_type)
  use jcup_interface, only : jcup_set_local_coef
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: send_comp
  integer, intent(IN) :: mapping_tag
  real(kind=8), intent(IN)    :: global_coef(:)
  real(kind=8), intent(INOUT) :: local_coef(:)
  integer, intent(IN) :: coef_type

  if (is_coupled) then
    call jcup_set_local_coef(my_comp, send_comp, mapping_tag, global_coef, local_coef, coef_type)
  end if

end subroutine h3ou_set_local_coef

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_coef(my_comp, recv_comp, coef)
  use jcup_interface, only : jcup_send_coef
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: recv_comp
  real(kind=8), intent(IN)     :: coef(:)

  if (is_coupled) call jcup_send_coef(my_comp, recv_comp, coef)

end subroutine h3ou_send_coef

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_coef(my_comp, send_comp, coef, coef_type)
  use jcup_interface, only : jcup_recv_coef
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: send_comp
  real(kind=8), intent(INOUT)  :: coef(:)
  integer, intent(IN) :: coef_type

  if (is_coupled) call jcup_recv_coef(my_comp, send_comp, 1, coef, coef_type)

end subroutine h3ou_recv_coef


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_get_local_operation_index(recv_model_name, send_model_name, grid_tag, &
                                     num_of_operation, operation_index, &
                                     send_data_index, recv_data_index, &
                                     send_coef_index, recv_coef_index)
  use jcup_interface, only : jcup_get_local_operation_index
  implicit none
  character(len=*), intent(IN) :: recv_model_name, send_model_name
  integer, optional, intent(IN) :: grid_tag
  integer, intent(INOUT) :: num_of_operation
  integer, pointer :: operation_index(:)
  integer, pointer :: send_data_index(:)
  integer, pointer :: recv_data_index(:)
  integer, pointer :: send_coef_index(:)
  integer, pointer :: recv_coef_index(:)

  if (is_coupled) then
     call jcup_get_local_operation_index(recv_model_name, send_model_name, grid_tag, &
       num_of_operation, operation_index, &
       send_data_index, recv_data_index, &
       send_coef_index, recv_coef_index)
  end if
  
end subroutine h3ou_get_local_operation_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_get_num_of_send_grid(recv_model_name, send_model_name, grid_tag, num_of_grid)
  use jcup_interface, only : jcup_get_num_of_send_grid
  implicit none
  character(len=*), intent(IN) :: recv_model_name, send_model_name
  integer, optional, intent(IN) :: grid_tag
  integer, intent(INOUT) :: num_of_grid

  if (is_coupled) then
     call jcup_get_num_of_send_grid(recv_model_name, send_model_name, grid_tag, num_of_grid)
  end if
  
end subroutine h3ou_get_num_of_send_grid

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_get_num_of_recv_grid(recv_model_name, send_model_name, grid_tag, num_of_grid)
  use jcup_interface, only : jcup_get_num_of_recv_grid
  implicit none
  character(len=*), intent(IN) :: recv_model_name, send_model_name
  integer, optional, intent(IN) :: grid_tag
  integer, intent(INOUT) :: num_of_grid

  if (is_coupled) then
     call jcup_get_num_of_recv_grid(recv_model_name, send_model_name, grid_tag, num_of_grid)
  end if
  
end subroutine h3ou_get_num_of_recv_grid

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_init_time(time_array)
  use mod_utils, only : put_log, close_log_file
  use jcup_interface, only : jcup_init_time
  use mod_utils, only : put_log
  use jcup_mpi_lib, only : jml_finalize
  implicit none
  integer, intent(IN) :: time_array(6)
  character(len=STR_LONG) :: log_str
  
  if (.not.is_coupled) return
  
  call put_log("")
  call put_log(">>>>>>>>>>>>>>> h3ou_init_time IN")
  call put_log("--------------- initial time information")
  write(log_str, '(A,I04,"/"I02,"/",I02,"-",I02,":",I02,":",I02)') " initial_time = ",time_array(1), time_array(2), time_array(3), &
                                                                                      time_array(4), time_array(5), time_array(6)
  call put_log(trim(log_str))

  if (stop_step == 7) then
    call put_log("[stop_step = 7], h3ou_init_time:IN   OK")
    call close_log_file()
    call jml_finalize(.true.)
    write(0, *) "77777 --- Coupler stopped at the beginning of h3ou_init_time --- 77777"
    stop 
  end if

  call jcup_init_time(time_array)

  call put_log("<<<<<<<<<<<<<<< h3ou_init_time OUT")

  if (stop_step == 8) then
    call put_log("[stop_step = 8], h3ou_init_time:OUT  OK. Breakpoint check run completed.")
    call close_log_file()
    call jml_finalize(.true.)
    write(0, *) "88888 --- Coupler stopped at the end of h3ou_init_time --- 88888"
    stop 
  end if

end subroutine h3ou_init_time

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_set_time(component_name, time_array, delta_t, is_exchange)
  use mod_utils, only : put_log, close_log_file
  use jcup_interface, only : jcup_set_time, jcup_log
  use palmtime, only : palm_TimeStart, palm_TimeEnd
  use mod_ensemble_base, only : NO_ENSEMBLE, get_ensemble_type, get_my_name
  use mod_utils, only : put_log
  use jcup_mpi_lib, only : jml_finalize
  implicit none
  character(len=*), intent(IN) :: component_name
  integer, intent(IN)     :: time_array(:) ! 2014/07/03 [MOD]
  integer, intent(IN)     :: delta_t
  logical, optional       :: is_exchange
  character(len=STR_LONG) :: log_str
  integer :: ierr


  call palm_TimeStart("h3ou_set_time")
  
  call put_log("",1)
  call put_log(">>>>>>>>>>>>>>> h3ou_set_time IN",1)
  call put_log("--------------- current time information",1)
  write(log_str, '(A,I04,"/"I02,"/",I02,"-",I02,":",I02,":",I02, A, I0)') &
                 " current_time = ",time_array(1), time_array(2), time_array(3), &
                 time_array(4), time_array(5), time_array(6), &
                 ", delta_t = ", delta_t
  call put_log(trim(log_str),1)
 
  if (stop_step == 9) then
    call put_log("[stop_step = 9], h3ou_set_time:IN   OK")
    call close_log_file()
    call jml_finalize(.true.)
    write(0, *) "99999 --- Coupler stopped at the beginning of h3ou_set_time --- 99999"
    stop 
  end if

  if (is_coupled) then
    if (get_ensemble_type() == NO_ENSEMBLE) then
       call jcup_set_time(component_name, time_array, delta_t, is_exchange)
    else
       call jcup_set_time(trim(get_my_name()), time_array, delta_t, is_exchange)
    end if
  else
     write(log_str, '(A,I0,"-",I0,"-",I0," ",I0,":",I0,":",I0, A, I0)') &
          "current_time : ",time_array(1), time_array(2), &
           time_array(3), time_array(4), time_array(5), time_array(6), &
                                                     ",  delat T = ", delta_t
    call jcup_log("h3ou_set_time", trim(log_str), 2)
  end if

  call palm_TimeEnd("h3ou_set_time")

  call put_log("<<<<<<<<<<<<<<< h3ou_set_time OUT",1)

  if (stop_step == 10) then
    call put_log("[stop_step = 10], h3ou_set_time:OUT  OK. Breakpoint check run completed.")
    call close_log_file()
    call jml_finalize(.true.)
    write(0, *) "101010 --- Coupler stopped at the end of h3ou_set_time --- 101010"
    stop 
  end if

end subroutine h3ou_set_time

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_inc_time(time_array, delta_t)
  use jcup_interface, only : jcup_inc_calendar
  implicit none
  integer, intent(INOUT) :: time_array(:)
  integer, intent(IN)    :: delta_t
  
  if (is_coupled) then
    call jcup_inc_calendar(time_array, delta_t)
  end if
 
end subroutine h3ou_inc_time

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_put_data_1d(data_name, data, data_vector)
  use jcup_interface, only : jcup_put_data, jcup_log
  use mod_utils, only : put_error, put_log
  use mod_comp, only : comp_type, get_comp_ptr, varp_type, get_varp_ptr
  use mod_ensemble_base, only : get_ensemble_type, get_ensemble_mean
  use palmtime, only : palm_TimeStart, palm_TimeEnd
  implicit none
  character(len=*), intent(IN) :: data_name
  real(kind=8), intent(INOUT) :: data(:)
  real(kind=8), optional, intent(IN) :: data_vector(:)
  type(comp_type), pointer :: comp_ptr
  type(varp_type), pointer :: varp_ptr
  logical                  :: is_get_ok
  character(len=STR_LONG)  :: log_str
  
  call put_log("",2)
  call put_log(">>>>>>>>>>>>>>> h3ou_put_data IN",2)
  call put_log("--------------- put data information",2)
  write(log_str, '(A,F20.8,A,F20.8)') " put_data_name = "//trim(data_name)//", min = ", minval(data), ", max = ", maxval(data)
  call put_log(trim(log_str), 2)
  
  call palm_TimeStart("put_data_1d")
  
  comp_ptr => get_comp_ptr(trim(my_name))
  varp_ptr => get_varp_ptr(comp_ptr%pdata_ptr, trim(data_name))

  if (.not.associated(varp_ptr)) then
     call put_error("h3ou_put_data_1d", "Data "//trim(data_name)//" is not defined in config file")
  end if

  if (get_ensemble_type() == NO_ENSEMBLE) then
     call jcup_put_data(varp_ptr%jcup_varp_ptr, data, data_vector)
     call palm_TimeEnd("put_data_1d")
     return
  end if
  
  if (is_coupled) then
     if (get_ensemble_type() == ONE_TO_ONE) then
      call jcup_put_data(varp_ptr%jcup_varp_ptr, data, data_vector)
    else ! root of MANY_TO_ONE
       call get_ensemble_mean(data, is_get_ok)
       call jcup_log("h3ou_put_data", "cal ensemble mean", 2)
       call jcup_put_data(varp_ptr%jcup_varp_ptr, data, data_vector)     
    end if
  else ! no root of MANY_TO_ONE
    call jcup_log("h3ou_put_data", "put and cal ensemble mean", 2)
    call get_ensemble_mean(data, is_get_ok)
  end if
    
  call palm_TimeEnd("put_data_1d")

  call put_log("<<<<<<<<<<<<<<< h3ou_put_data OUT",2)

end subroutine h3ou_put_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_put_data_25d(data_name, data, data_vector)
  use jcup_interface, only : jcup_put_data, jcup_log
  use mod_utils, only : put_error, put_log
  use mod_comp, only : comp_type, get_comp_ptr, varp_type, get_varp_ptr
  use mod_ensemble_base, only : get_ensemble_type, get_ensemble_mean
  use palmtime, only : palm_TimeStart, palm_TimeEnd
  implicit none
  character(len=*), intent(IN) :: data_name
  real(kind=8), intent(INOUT) :: data(:,:)
  real(kind=8), optional, intent(IN) :: data_vector(:)
  type(comp_type), pointer :: comp_ptr
  type(varp_type), pointer :: varp_ptr
  logical                  :: is_get_ok
  character(len=STR_LONG)  :: log_str
  
  call put_log("",2)
  call put_log(">>>>>>>>>>>>>>> h3ou_put_data IN",2)
  call put_log("--------------- put data information",2)
  write(log_str, '(A,F20.8,A,F20.8)') " put_data_name = "//trim(data_name)//", min = ", minval(data), ", max = ", maxval(data)
  call put_log(trim(log_str), 2)
  
  call palm_TimeStart("put_data_2d")

  comp_ptr => get_comp_ptr(trim(my_name))
  varp_ptr => get_varp_ptr(comp_ptr%pdata_ptr, trim(data_name))

  if (.not.associated(varp_ptr)) then
     call put_error("h3ou_put_data_25d", "Data "//trim(data_name)//" is not defined in config file")
  end if

  if (get_ensemble_type() == NO_ENSEMBLE) then
     call jcup_put_data(varp_ptr%jcup_varp_ptr, data, data_vector)
     call palm_TimeEnd("put_data_2d")
     return
  end if

  if (is_coupled) then
    if (get_ensemble_type() == ONE_TO_ONE) then
       call jcup_put_data(varp_ptr%jcup_varp_ptr, data, data_vector)
    else ! root of MANY_TO_ONE
       call get_ensemble_mean(data, is_get_ok)
       call jcup_log("h3ou_put_data", "cal ensemble mean", 2)
       call jcup_put_data(varp_ptr%jcup_varp_ptr, data, data_vector)     
    end if
  else ! no root of MANY_TO_ONE
    call jcup_log("h3ou_put_data", "put and cal ensemble mean", 2)
    call get_ensemble_mean(data, is_get_ok)
  end if

  call palm_TimeEnd("put_data_2d")

  call put_log("<<<<<<<<<<<<<<< h3ou_put_data OUT",2)

end subroutine h3ou_put_data_25d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_get_data_1d(data_name, data, data_vector, is_recv_ok)
  use jcup_interface, only : jcup_get_data, jcup_log
  use mod_utils, only : put_error, put_log
  use mod_comp, only : comp_type, get_comp_ptr, varg_type, get_varg_ptr
  use mod_ensemble_base, only : get_ensemble_type, bcast_ensemble
  use palmtime, only : palm_TimeStart, palm_TimeEnd
  implicit none
  character(len=*), intent(IN) :: data_name
  real(kind=8), intent(INOUT) :: data(:)
  real(kind=8), intent(OUT), optional :: data_vector(:)
  type(comp_type), pointer :: comp_ptr
  type(varg_type), pointer :: varg_ptr
  logical, intent(OUT)     :: is_recv_ok
  integer                  :: recv_flag
  character(len=STR_LONG)  :: log_str

  call put_log("",2)
  call put_log(">>>>>>>>>>>>>>> h3ou_get_data IN",2)

  call palm_TimeStart("get_data_1d")

  comp_ptr => get_comp_ptr(trim(my_name))
  varg_ptr => get_varg_ptr(comp_ptr%gdata_ptr, trim(data_name))

  
  if (.not.associated(varg_ptr)) then
     call put_error("h3ou_get_data_1d", "Data "//trim(data_name)//" is not defined in config file")
  end if
  
  if (get_ensemble_type() == NO_ENSEMBLE) then
     call jcup_get_data(varg_ptr%jcup_varg_ptr, data, data_vector, is_recv_ok = is_recv_ok)
     goto 8000
  end if

  if (is_coupled) then
     if (get_ensemble_type() == ONE_TO_ONE) then
        call jcup_get_data(varg_ptr%jcup_varg_ptr, data, data_vector, is_recv_ok = is_recv_ok)
     else ! root of MANY_TO_ONE
        call jcup_get_data(varg_ptr%jcup_varg_ptr, data, data_vector, is_recv_ok = is_recv_ok)
        if (is_recv_ok) then
           Recv_flag = 1
        else
           recv_flag = 0
        end if
        call bcast_ensemble(recv_flag)
        if (is_recv_ok) then 
           call bcast_ensemble(data)
           call jcup_log("h3ou_get_data", "bcast data for ensemble", 2)
         end if
     end if
  else ! no root of MANY_TO_ONE
     call bcast_ensemble(recv_flag)
     if (recv_flag == 1) then
        is_recv_ok = .true.
     else
        is_recv_ok = .false.
     end if
     if (is_recv_ok) then 
        call bcast_ensemble(data)
        call jcup_log("h3ou_get_data", "get ensemble data", 2)
     end if
  end if
  
  8000 continue

  call palm_TimeEnd("get_data_1d")
  
  call put_log("--------------- get data information",2)
  if (is_recv_ok) then
    write(log_str, '(A,F20.8,A,F20.8)') " get_data_name = "//trim(data_name)//", min = ", minval(data), ", max = ", maxval(data)
    call put_log(trim(log_str), 2)
  else
    write(log_str, '(A)') " get_data_name = "//trim(data_name)//", not get step"
    call put_log(trim(log_str), 2)
  end if

  call put_log("<<<<<<<<<<<<<<< h3ou_get_data OUT",2)

end subroutine h3ou_get_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_get_data_25d(data_name, data, data_vector, is_recv_ok)
  use jcup_interface, only : jcup_get_data, jcup_log
  use mod_utils, only : put_error, put_log
  use mod_comp, only : comp_type, get_comp_ptr, varg_type, get_varg_ptr
  use mod_ensemble_base, only : get_ensemble_type, bcast_ensemble
  use palmtime, only : palm_TimeStart, palm_TimeEnd
  implicit none
  character(len=*), intent(IN) :: data_name
  real(kind=8), intent(INOUT) :: data(:,:)
  real(kind=8), intent(OUT), optional :: data_vector(:)
  type(comp_type), pointer :: comp_ptr
  type(varg_type), pointer :: varg_ptr
  logical, intent(OUT)     :: is_recv_ok
  integer                  :: recv_flag
  character(len=STR_LONG)  :: log_str
  
  call put_log("",2)
  call put_log(">>>>>>>>>>>>>>> h3ou_get_data IN",2)

  call palm_TimeStart("get_data_2d")

  comp_ptr => get_comp_ptr(trim(my_name))
  varg_ptr => get_varg_ptr(comp_ptr%gdata_ptr, trim(data_name))

  if (.not.associated(varg_ptr)) then
     call put_error("h3ou_get_data_25d", "Data "//trim(data_name)//" is not defined in config file")
  end if

  if (get_ensemble_type() == NO_ENSEMBLE) then
     call jcup_get_data(varg_ptr%jcup_varg_ptr, data, data_vector, is_recv_ok = is_recv_ok)
     goto 8000
  end if

  if (is_coupled) then
     if (get_ensemble_type() == ONE_TO_ONE) then
        call jcup_get_data(varg_ptr%jcup_varg_ptr, data, data_vector, is_recv_ok)
     else ! root of MANY_TO_ONE
        call jcup_get_data(varg_ptr%jcup_varg_ptr, data, data_vector, is_recv_ok)
        if (is_recv_ok) then
           recv_flag = 1
        else
           recv_flag = 0
        end if
        call bcast_ensemble(recv_flag)
        if (is_recv_ok) then
           call bcast_ensemble(data)
           call jcup_log("h3ou_get_data", "bcast data for ensemble", 2)
        end if
     end if
  else ! no root of MANY_TO_ONE
     call bcast_ensemble(recv_flag)
     if (recv_flag == 1) then
        is_recv_ok = .true.
     else
        is_recv_ok = .false.
     end if
     if (is_recv_ok) then 
        call bcast_ensemble(data)
        call jcup_log("h3ou_get_data", "get ensemble data", 2)
     end if
  end if
  
  8000 continue

  call palm_TimeEnd("get_data_2d")

  call put_log("--------------- get data information",2)
  if (is_recv_ok) then
    write(log_str, '(A,F20.8,A,F20.8)') " get_data_name = "//trim(data_name)//", min = ", minval(data), ", max = ", maxval(data)
    call put_log(trim(log_str), 2)
  else
    write(log_str, '(A)') " get_data_name = "//trim(data_name)//", not get step"
    call put_log(trim(log_str), 2)
  end if

  call put_log("<<<<<<<<<<<<<<< h3ou_get_data OUT",2)

end subroutine h3ou_get_data_25d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_data_immediately(send_comp, recv_comp, time_lag)
  use jcup_interface, only : jcup_send_data_immediately
  implicit none
  character(len=*), intent(IN) :: send_comp 
  character(len=*), intent(IN) :: recv_comp 
  integer, intent(IN) :: time_lag

  if (is_coupled) call jcup_send_data_immediately(send_comp, recv_comp, time_lag)

end subroutine h3ou_send_data_immediately

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_data_immediately(send_comp, recv_comp)
  use jcup_interface, only : jcup_recv_data_immediately
  implicit none
  character(len=*), intent(IN) :: send_comp 
  character(len=*), intent(IN) :: recv_comp 

  if (is_coupled) call jcup_recv_data_immediately(send_comp, recv_comp)

end subroutine h3ou_recv_data_immediately


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_array_str(my_comp, recv_comp, array)
  use jcup_interface, only : jcup_send_array
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: recv_comp
  character(len=*), intent(IN) :: array

  if (is_coupled) call jcup_send_array(my_comp, recv_comp, array)

end subroutine h3ou_send_array_str

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_array_int(my_comp, recv_comp, array)
  use jcup_interface, only : jcup_send_array
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: recv_comp
  integer, intent(IN) :: array(:)

  if (is_coupled) call jcup_send_array(my_comp, recv_comp, array)

end subroutine h3ou_send_array_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_array_real(my_comp, recv_comp, array)
  use jcup_interface, only : jcup_send_array
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: recv_comp
  real(kind=4), intent(IN) :: array(:)

  if (is_coupled) call jcup_send_array(my_comp, recv_comp, array)

end subroutine h3ou_send_array_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_array_dbl(my_comp, recv_comp, array)
  use jcup_interface, only : jcup_send_array
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: recv_comp
  real(kind=8), intent(IN) :: array(:)

  if (is_coupled) call jcup_send_array(my_comp, recv_comp, array)

end subroutine h3ou_send_array_dbl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_array_str(my_comp, send_comp, array)
  use jcup_interface, only : jcup_recv_array
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: send_comp
  character(len=*), intent(INOUT) :: array

  if (is_coupled) call jcup_recv_array(my_comp, send_comp, array)

end subroutine h3ou_recv_array_str

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_array_int(my_comp, send_comp, array)
  use jcup_interface, only : jcup_recv_array
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: send_comp
  integer, intent(INOUT) :: array(:)

  if (is_coupled) call jcup_recv_array(my_comp, send_comp, array)

end subroutine h3ou_recv_array_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_array_real(my_comp, send_comp, array)
  use jcup_interface, only : jcup_recv_array
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: send_comp
  real(kind=4), intent(INOUT) :: array(:)

  if (is_coupled) call jcup_recv_array(my_comp, send_comp, array)

end subroutine h3ou_recv_array_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_array_dbl(my_comp, send_comp, array)
  use jcup_interface, only : jcup_recv_array
  implicit none
  character(len=*), intent(IN) :: my_comp
  character(len=*), intent(IN) :: send_comp
  real(kind=8), intent(INOUT) :: array(:)

  if (is_coupled) call jcup_recv_array(my_comp, send_comp, array)

end subroutine h3ou_recv_array_dbl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_array_str_2(recv_comp, array)
  use jcup_interface, only : jcup_send_array
  implicit none
  character(len=*), intent(IN) :: recv_comp
  character(len=*), intent(IN) :: array

  if (is_coupled) call jcup_send_array(my_name, recv_comp, array)

end subroutine h3ou_send_array_str_2

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_array_int_2(recv_comp, array)
  use jcup_interface, only : jcup_send_array
  implicit none
  character(len=*), intent(IN) :: recv_comp
  integer, intent(IN) :: array(:)

  if (is_coupled) call jcup_send_array(my_name, recv_comp, array)

end subroutine h3ou_send_array_int_2

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_array_real_2(recv_comp, array)
  use jcup_interface, only : jcup_send_array
  implicit none
  character(len=*), intent(IN) :: recv_comp
  real(kind=4), intent(IN) :: array(:)

  if (is_coupled) call jcup_send_array(my_name, recv_comp, array)

end subroutine h3ou_send_array_real_2

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_array_dbl_2(recv_comp, array)
  use jcup_interface, only : jcup_send_array
  implicit none
  character(len=*), intent(IN) :: recv_comp
  real(kind=8), intent(IN) :: array(:)

  if (is_coupled) call jcup_send_array(my_name, recv_comp, array)

end subroutine h3ou_send_array_dbl_2

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_array_str_2(send_comp, array)
  use jcup_interface, only : jcup_recv_array
  implicit none
  character(len=*), intent(IN) :: send_comp
  character(len=*), intent(INOUT) :: array

  if (is_coupled) call jcup_recv_array(my_name, send_comp, array)

end subroutine h3ou_recv_array_str_2

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_array_int_2(send_comp, array)
  use jcup_interface, only : jcup_recv_array
  implicit none
  character(len=*), intent(IN) :: send_comp
  integer, intent(INOUT) :: array(:)

  if (is_coupled) call jcup_recv_array(my_name, send_comp, array)

end subroutine h3ou_recv_array_int_2

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_array_real_2(send_comp, array)
  use jcup_interface, only : jcup_recv_array
  implicit none
  character(len=*), intent(IN) :: send_comp
  real(kind=4), intent(INOUT) :: array(:)

  if (is_coupled) call jcup_recv_array(my_name, send_comp, array)

end subroutine h3ou_recv_array_real_2

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_array_dbl_2(send_comp, array)
  use jcup_interface, only : jcup_recv_array
  implicit none
  character(len=*), intent(IN) :: send_comp
  real(kind=8), intent(INOUT) :: array(:)

  if (is_coupled) call jcup_recv_array(my_name, send_comp, array)

end subroutine h3ou_recv_array_dbl_2

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+            Simple Mode APIs           +=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+


subroutine h3ou_end()
  implicit none

  call mpi_finalize(ierror)

end subroutine h3ou_end

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+            Global Routines            +=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_global_string(source_name, str_data)
  use jcup_interface, only : jcup_get_model_id
  use jcup_mpi_lib, only : jml_GetLeaderRank, jml_BcastGlobal
  implicit none
  character(len=*), intent(IN)    :: source_name
  character(len=*), intent(INOUT) :: str_data
  integer :: model_id, source_pe

  call jcup_get_model_id(trim(source_name), model_id)
  
  source_pe = jml_GetLeaderRank(model_id)

  call jml_BcastGlobal(str_data, source_pe)

end subroutine h3ou_bcast_global_string

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_global_int(source_name, data)
  use jcup_interface, only : jcup_get_model_id
  use jcup_mpi_lib, only : jml_GetLeaderRank, jml_BcastGlobal
  implicit none
  character(len=*), intent(IN)    :: source_name
  integer, intent(INOUT) :: data(:)
  integer :: model_id, source_pe

  call jcup_get_model_id(trim(source_name), model_id)
  
  source_pe = jml_GetLeaderRank(model_id)

  call jml_BcastGlobal(data, 1, size(data), source_pe)

end subroutine h3ou_bcast_global_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_global_real(source_name, data)
  use jcup_interface, only : jcup_get_model_id
  use jcup_mpi_lib, only : jml_GetLeaderRank, jml_BcastGlobal
  implicit none
  character(len=*), intent(IN)    :: source_name
  real(kind=4), intent(INOUT) :: data(:)
  integer :: model_id, source_pe

  call jcup_get_model_id(trim(source_name), model_id)
  
  source_pe = jml_GetLeaderRank(model_id)

  call jml_BcastGlobal(data, 1, size(data), source_pe)

end subroutine h3ou_bcast_global_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_global_double(source_name, data)
  use jcup_interface, only : jcup_get_model_id
  use jcup_mpi_lib, only : jml_GetLeaderRank, jml_BcastGlobal
  implicit none
  character(len=*), intent(IN)    :: source_name
  real(kind=8), intent(INOUT) :: data(:)
  integer :: model_id, source_pe

  call jcup_get_model_id(trim(source_name), model_id)
  
  source_pe = jml_GetLeaderRank(model_id)

  call jml_BcastGlobal(data, 1, size(data), source_pe)

end subroutine h3ou_bcast_global_double

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+           Inter Model Routines        +=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_model_int(source_name, target_name, data)
  use jcup_mpi_lib, only : jml_isLocalLeader, jml_BcastLocal, jml_SendLeader, jml_RecvLeader
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: source_name
  character(len=*), intent(IN) :: target_name
  integer, intent(INOUT)       :: data(:)
  integer :: send_id, recv_id
  
  if (trim(source_name) == trim(my_name)) then
     if (jml_isLocalLeader(my_id)) then
        call jcup_get_model_id(trim(target_name), recv_id)
        call jml_SendLeader(data, 1, size(data), recv_id-1)
     end if
     call jml_BcastLocal(my_id, data, 1, size(data))
  end if
     
  if (trim(target_name) == trim(my_name)) then
     if (jml_isLocalLeader(my_id)) then
        call jcup_get_model_id(trim(source_name), send_id)
        call jml_RecvLeader(data, 1, size(data), send_id-1)
     end if
     call jml_BcastLocal(my_id, data, 1, size(data))
  end if
  
end subroutine h3ou_bcast_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_model_real(source_name, target_name, data)
  use jcup_mpi_lib, only : jml_isLocalLeader, jml_BcastLocal, jml_SendLeader, jml_RecvLeader
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: source_name
  character(len=*), intent(IN) :: target_name
  real(kind=4), intent(INOUT)  :: data(:)
  integer :: send_id, recv_id
  
  if (trim(source_name) == trim(my_name)) then
     if (jml_isLocalLeader(my_id)) then
        call jcup_get_model_id(trim(target_name), recv_id)
        call jml_SendLeader(data, 1, size(data), recv_id-1)
     end if
     call jml_BcastLocal(my_id, data, 1, size(data))
  end if
     
  if (trim(target_name) == trim(my_name)) then
     if (jml_isLocalLeader(my_id)) then
        call jcup_get_model_id(trim(source_name), send_id)
        call jml_RecvLeader(data, 1, size(data), send_id-1)
     end if
     call jml_BcastLocal(my_id, data, 1, size(data))
  end if
  
end subroutine h3ou_bcast_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_model_double(source_name, target_name, data)
  use jcup_mpi_lib, only : jml_isLocalLeader, jml_BcastLocal, jml_SendLeader, jml_RecvLeader
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: source_name
  character(len=*), intent(IN) :: target_name
  real(kind=8), intent(INOUT)  :: data(:)
  integer :: send_id, recv_id
  
  if (trim(source_name) == trim(my_name)) then
     if (jml_isLocalLeader(my_id)) then
        call jcup_get_model_id(trim(target_name), recv_id)
        call jml_SendLeader(data, 1, size(data), recv_id-1)
     end if
     call jml_BcastLocal(my_id, data, 1, size(data))
  end if
     
  if (trim(target_name) == trim(my_name)) then
     if (jml_isLocalLeader(my_id)) then
        call jcup_get_model_id(trim(source_name), send_id)
        call jml_RecvLeader(data, 1, size(data), send_id-1)
     end if
     call jml_BcastLocal(my_id, data, 1, size(data))
  end if
  
end subroutine h3ou_bcast_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_model_int(target_name, target_pe, data)
  use jcup_mpi_lib, only : jml_ISendModel, jml_send_waitall
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN)          :: target_pe
  integer, target, intent(IN)  :: data(:)
  integer :: target_id
  integer, pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(target_name), target_id)

  call jml_ISendModel(my_id, data_ptr, 1, size(data), target_id, target_pe)
  call jml_send_waitall()
  
end subroutine h3ou_send_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_model_real(target_name, target_pe, data)
  use jcup_mpi_lib, only : jml_ISendModel, jml_send_waitall
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN)          :: target_pe
  real(kind=4), target, intent(IN)  :: data(:)
  integer :: target_id
  real(kind=4), pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(target_name), target_id)

  call jml_ISendModel(my_id, data_ptr, 1, size(data), target_id, target_pe)
  call jml_send_waitall()
  
end subroutine h3ou_send_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_model_double(target_name, target_pe, data)
  use jcup_mpi_lib, only : jml_ISendModel, jml_send_waitall
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN)          :: target_pe
  real(kind=8), target, intent(IN)  :: data(:)
  integer :: target_id
  real(kind=8), pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(target_name), target_id)

  call jml_ISendModel(my_id, data_ptr, 1, size(data), target_id, target_pe)
  call jml_send_waitall()
  
end subroutine h3ou_send_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_model_int(source_name, source_pe, data)
  use jcup_mpi_lib, only : jml_IRecvModel, jml_recv_waitall
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN)   :: source_name
  integer, intent(IN)            :: source_pe
  integer, target, intent(INOUT) :: data(:)
  integer :: source_id
  integer, pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(source_name), source_id)

  call jml_IRecvModel(my_id, data_ptr, 1, size(data), source_id, source_pe)
  call jml_recv_waitall()
  
end subroutine h3ou_recv_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_model_real(source_name, source_pe, data)
  use jcup_mpi_lib, only : jml_IRecvModel, jml_recv_waitall
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN)          :: source_name
  integer, intent(IN)                  :: source_pe
  real(kind=4), target, intent(INOUt)  :: data(:)
  integer :: source_id
  real(kind=4), pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(source_name), source_id)

  call jml_IRecvModel(my_id, data_ptr, 1, size(data), source_id, source_pe)
  call jml_recv_waitall()
  
end subroutine h3ou_recv_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_model_double(source_name, source_pe, data)
  use jcup_mpi_lib, only : jml_IRecvModel, jml_recv_waitall
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN)          :: source_name
  integer, intent(IN)                  :: source_pe
  real(kind=8), target, intent(INOUt)  :: data(:)
  integer :: source_id
  real(kind=8), pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(source_name), source_id)

  call jml_IRecvModel(my_id, data_ptr, 1, size(data), source_id, source_pe)
  call jml_recv_waitall()
  
end subroutine h3ou_recv_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_int_scalar(target_name, val)
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN) :: val
  integer :: send_array(1)

  send_array(1) = val

  call h3ou_send_array(my_name, target_name, send_array)
  
end subroutine h3ou_send_int_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_real_scalar(target_name, val)
  implicit none
  character(len=*), intent(IN) :: target_name
  real(kind=4), intent(IN) :: val
  real(kind=4) :: send_array(1)

  send_array(1) = val

  call h3ou_send_array(my_name, target_name, send_array)
  
end subroutine h3ou_send_real_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_double_scalar(target_name, val)
  implicit none
  character(len=*), intent(IN) :: target_name
  real(kind=8), intent(IN) :: val
  real(kind=8) :: send_array(1)

  send_array(1) = val

  call h3ou_send_array(my_name, target_name, send_array)
  
end subroutine h3ou_send_double_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_int_scalar(source_name, val)
  implicit none
  character(len=*), intent(IN) :: source_name
  integer, intent(INOUT) :: val
  integer :: recv_array(1)


  call h3ou_recv_array(my_name, source_name, recv_array)
  
  val = recv_array(1)
  
end subroutine h3ou_recv_int_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_real_scalar(source_name, val)
  implicit none
  character(len=*), intent(IN) :: source_name
  real(kind=4), intent(INOUT) :: val
  real(kind=4) :: recv_array(1)


  call h3ou_recv_array(my_name, source_name, recv_array)
  
  val = recv_array(1)
  
end subroutine h3ou_recv_real_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_double_scalar(source_name, val)
  implicit none
  character(len=*), intent(IN) :: source_name
  real(kind=8), intent(INOUT) :: val
  real(kind=8) :: recv_array(1)

  call h3ou_recv_array(my_name, source_name, recv_array)
  
  val = recv_array(1)
  
end subroutine h3ou_recv_double_scalar

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_isend_model_int(target_name, target_pe, data)
  use jcup_mpi_lib, only : jml_ISendModel
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN)          :: target_pe
  integer, target, intent(IN)  :: data(:)
  integer :: target_id
  integer, pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(target_name), target_id)

  call jml_ISendModel(my_id, data_ptr, 1, size(data), target_id, target_pe)
  
end subroutine h3ou_isend_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_isend_model_real(target_name, target_pe, data)
  use jcup_mpi_lib, only : jml_ISendModel
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN)          :: target_pe
  real(kind=4), target, intent(IN)  :: data(:)
  integer :: target_id
  real(kind=4), pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(target_name), target_id)

  call jml_ISendModel(my_id, data_ptr, 1, size(data), target_id, target_pe)
  
end subroutine h3ou_isend_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_isend_model_double(target_name, target_pe, data)
  use jcup_mpi_lib, only : jml_ISendModel
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN) :: target_name
  integer, intent(IN)          :: target_pe
  real(kind=8), target, intent(IN)  :: data(:)
  integer :: target_id
  real(kind=8), pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(target_name), target_id)

  call jml_ISendModel(my_id, data_ptr, 1, size(data), target_id, target_pe)
  
end subroutine h3ou_isend_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_irecv_model_int(source_name, source_pe, data)
  use jcup_mpi_lib, only : jml_IRecvModel
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN)   :: source_name
  integer, intent(IN)            :: source_pe
  integer, target, intent(INOUT) :: data(:)
  integer :: source_id
  integer, pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(source_name), source_id)

  call jml_IRecvModel(my_id, data_ptr, 1, size(data), source_id, source_pe)
  
end subroutine h3ou_irecv_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_irecv_model_real(source_name, source_pe, data)
  use jcup_mpi_lib, only : jml_IRecvModel
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN)          :: source_name
  integer, intent(IN)                  :: source_pe
  real(kind=4), target, intent(INOUt)  :: data(:)
  integer :: source_id
  real(kind=4), pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(source_name), source_id)

  call jml_IRecvModel(my_id, data_ptr, 1, size(data), source_id, source_pe)
  
end subroutine h3ou_irecv_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_irecv_model_double(source_name, source_pe, data)
  use jcup_mpi_lib, only : jml_IRecvModel
  use jcup_interface, only : jcup_get_model_id
  implicit none
  character(len=*), intent(IN)          :: source_name
  integer, intent(IN)                  :: source_pe
  real(kind=8), target, intent(INOUt)  :: data(:)
  integer :: source_id
  real(kind=8), pointer :: data_ptr

  data_ptr => data(1)
  
  call jcup_get_model_id(trim(source_name), source_id)

  call jml_IRecvModel(my_id, data_ptr, 1, size(data), source_id, source_pe)
  
end subroutine h3ou_irecv_model_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_isend_waitall()
  use jcup_mpi_lib, only : jml_send_waitall
  implicit none
  
  call jml_send_waitall()

end subroutine h3ou_isend_waitall


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_irecv_waitall()
  use jcup_mpi_lib, only : jml_recv_waitall
  implicit none
  
  call jml_recv_waitall()

end subroutine h3ou_irecv_waitall


!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+              Local Routines           +=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_local_int(source_pe, data)
  use jcup_mpi_lib, only : jml_BcastLocal
  implicit none
  integer, intent(IN)    :: source_pe
  integer, intent(INOUT) :: data(:)

  call jml_BcastLocal(my_id, data, 1, size(data), source_pe) 
  
end subroutine h3ou_bcast_local_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_local_real(source_pe, data)
  use jcup_mpi_lib, only : jml_BcastLocal
  implicit none
  integer, intent(IN)         :: source_pe
  real(kind=4), intent(INOUT) :: data(:)

  call jml_BcastLocal(my_id, data, 1, size(data), source_pe) 
  
end subroutine h3ou_bcast_local_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_bcast_local_double(source_pe, data)
  use jcup_mpi_lib, only : jml_BcastLocal
  implicit none
  integer, intent(IN)         :: source_pe
  real(kind=8), intent(INOUT) :: data(:)

  call jml_BcastLocal(my_id, data, 1, size(data), source_pe) 
  
end subroutine h3ou_bcast_local_double

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_local_int(target_pe, data)
  use jcup_mpi_lib, only : jml_SendLocal
  implicit none
  integer, intent(IN) :: target_pe
  integer, intent(IN) :: data(:)

  call jml_SendLocal(my_id, data, 1, size(data), target_pe)

end subroutine h3ou_send_local_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_local_real(target_pe, data)
  use jcup_mpi_lib, only : jml_SendLocal
  implicit none
  integer, intent(IN) :: target_pe
  real(kind=4), intent(IN) :: data(:)

  call jml_SendLocal(my_id, data, 1, size(data), target_pe)

end subroutine h3ou_send_local_real
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_send_local_double(target_pe, data)
  use jcup_mpi_lib, only : jml_SendLocal
  implicit none
  integer, intent(IN) :: target_pe
  real(kind=8), intent(IN) :: data(:)

  call jml_SendLocal(my_id, data, 1, size(data), target_pe)

end subroutine h3ou_send_local_double

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_local_int(source_pe, data)
  use jcup_mpi_lib, only : jml_RecvLocal
  implicit none
  integer, intent(IN) :: source_pe
  integer, intent(INOUT) :: data(:)

  call jml_RecvLocal(my_id, data, 1, size(data), source_pe)

end subroutine h3ou_recv_local_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_local_real(source_pe, data)
  use jcup_mpi_lib, only : jml_RecvLocal
  implicit none
  integer, intent(IN) :: source_pe
  real(kind=4), intent(INOUT) :: data(:)

  call jml_RecvLocal(my_id, data, 1, size(data), source_pe)

end subroutine h3ou_recv_local_real
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine h3ou_recv_local_double(source_pe, data)
  use jcup_mpi_lib, only : jml_RecvLocal
  implicit none
  integer, intent(IN) :: source_pe
  real(kind=8), intent(INOUT) :: data(:)

  call jml_RecvLocal(my_id, data, 1, size(data), source_pe)

end subroutine h3ou_recv_local_double
  
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

end module h3ou_api

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
! default interpolation 
subroutine interpolate_data(recv_model_name, send_model_name, mapping_tag, sn1, sn2, send_data, &
                            rn1, rn2, recv_data, num_of_data, tn, exchange_tag)
  use h3ou_api, only : h3ou_interpolation
  implicit none
  character(len=*), intent(IN) :: recv_model_name, send_model_name
  integer, intent(IN) :: mapping_tag
  integer, intent(IN) :: sn1, sn2
  real(kind=8), intent(IN) :: send_data(sn1,sn2)
  integer, intent(IN) :: rn1, rn2
  real(kind=8), intent(INOUT) :: recv_data(rn1,rn2)
  integer, intent(IN) :: num_of_data
  integer, intent(IN) :: tn
  integer, intent(IN) :: exchange_tag(tn)

  call h3ou_interpolation(recv_model_name, send_model_name, mapping_tag,&
       send_data, recv_data, num_of_data, exchange_tag(1))
  
end subroutine interpolate_data

