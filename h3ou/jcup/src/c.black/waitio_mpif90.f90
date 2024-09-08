! -*- f90 -*-
module waitio_mpif90
  implicit none
  include 'waitio_mpif.h'
  private
  !     -*- f90 -*-
  !     type definition

  ! /* -*- fortran -*- */

  !integer, public, parameter :: WAITIO_MPI_DEFTAG=524288
  !integer, public, parameter :: WAITIO_MPI_INTEGER=1+524288
  !integer, public, parameter :: WAITIO_MPI_INT=1+524288
  !integer, public, parameter :: WAITIO_MPI_LONG=2+524288
  !integer, public, parameter :: WAITIO_MPI_CHAR=3+524288
  !integer, public, parameter :: WAITIO_MPI_FLOAT=4+524288
  !integer, public, parameter :: WAITIO_MPI_DOUBLE=5+524288
  !integer, public, parameter :: WAITIO_MPI_DOUBLE_PRECISION=5+524288
  !integer, public, parameter :: WAITIO_MPI_SUM=1+524288
  !integer, public, parameter :: WAITIO_MPI_MAX=2+524288
  !integer, public, parameter :: WAITIO_MPI_MIN=3+524288
  !integer, public, parameter :: WAITIO_STATUS_SIZE=4
  !integer, public, parameter :: WAITIO_REQUEST_SIZE=22

  public :: WAITIO_MPI_Comm
  public :: WAITIO_MPI_Datatype
  public :: WAITIO_MPI_Op
  public :: WAITIO_MPI_Request
 
  public :: WMPI_init                    ! subroutine (ierror)
  public :: WMPI_Create_universe         ! subroutine (comm, ierror)
 
 
  type WAITIO_MPI_Datatype
     integer :: dtype
  end type WAITIO_MPI_Datatype

  type WAITIO_MPI_Op
     integer :: dop
  end type WAITIO_MPI_Op

  type WAITIO_MPI_Comm
     integer, pointer :: dpnt
  end type WAITIO_MPI_Comm

  type WAITIO_MPI_Request
     integer, dimension(*) :: dreq(22)
  end type WAITIO_MPI_Request

contains
  
  subroutine WMPI_Init (ierr)
    implicit none
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call WAITIO_INIT (1000, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WMPI_Init

  ! int waitio_create_universe (WAITIO_MPI_Comm *commp) ;

  subroutine WMPI_Create_universe ( comm, ierr)
    implicit none
    TYPE(WAITIO_MPI_Comm), INTENT(OUT) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call WAITIO_CREATE_UNIVERSE  ( comm, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WMPI_Create_Universe


  subroutine WAITIO_MPI_Isend(buf,count,datatype,dest,tag,comm,request,ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: buf
    INTEGER, INTENT(IN) :: count, dest, tag
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: datatype
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    TYPE(WAITIO_MPI_Request), INTENT(OUT) :: request
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_isend_f  ( buf, count, datatype, dest, tag, comm, request, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Isend

  ! int waitio_mpi_irecv ( void  *buf,  int count,  WAITIO_MPI_Datatype datatype,  int source,  int tag,
  !                        WAITIO_MPI_Comm comm,  WAITIO_MPI_Request  *request);
  subroutine WAITIO_MPI_Irecv ( buf, count, datatype, source, tag, comm, request, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: buf
    INTEGER, INTENT(IN) :: count, source, tag
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: datatype
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    TYPE(WAITIO_MPI_Request), INTENT(OUT) :: request
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_irecv_f  ( buf, count, datatype, source, tag, comm, request, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Irecv

  !int waitio_mpi_reduce ( const void  *sendbuf,  void  *recvbuf,  int count,  WAITIO_MPI_Datatype datatype,
  !                        WAITIO_MPI_Op op,  int root,  WAITIO_MPI_Comm comm) ;

  subroutine WAITIO_MPI_Reduce ( sendbuf, recvbuf, count, datatype, op, root, comm, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: sendbuf, recvbuf
    INTEGER, INTENT(IN) :: count, root
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: datatype
    TYPE(WAITIO_MPI_Op), INTENT(IN) :: op
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_reduce_f   ( sendbuf, recvbuf, count, datatype, op, root, comm, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Reduce

  !int waitio_mpi_bcast ( void  *buffer,  int count,  WAITIO_MPI_Datatype datatype,  int root,  
  !                       WAITIO_MPI_Comm comm) ;

  subroutine WAITIO_MPI_Bcast ( buffer, count, datatype, root, comm, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: buffer
    INTEGER, INTENT(IN) :: count, root
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: datatype
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_bcast_f   ( buffer, count, datatype, root, comm, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Bcast


  ! int waitio_mpi_allreduce ( const void  *sendbuf,  void  *recvbuf,  int count,  WAITIO_MPI_Datatype datatype,
  !                            WAITIO_MPI_Op op,  WAITIO_MPI_Comm comm) ;

  subroutine WAITIO_MPI_Allreduce ( sendbuf, recvbuf, count, datatype, op, comm, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: sendbuf, recvbuf
    INTEGER, INTENT(IN) :: count
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: datatype
    TYPE(WAITIO_MPI_Op), INTENT(IN) :: op
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_allreduce_f   ( sendbuf, recvbuf, count, datatype, op, comm, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Allreduce

  ! int waitio_mpi_waitall (int count,  WAITIO_MPI_Request *array_of_requests,  int  *array_of_statuses);

  subroutine WAITIO_MPI_Waitall (count, array_of_requests,  array_of_statuses, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: array_of_requests,  array_of_statuses
    INTEGER, INTENT(IN) :: count
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_waitall_f   (count, array_of_requests,  array_of_statuses, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Waitall


  !int waitio_mpi_gather ( const void  *sendbuf,  void  *recvbuf,  int count,  WAITIO_MPI_Datatype datatype,
  !                              int root,  WAITIO_MPI_Comm comm) ;

  subroutine WAITIO_MPI_Gather ( sendbuf, recvbuf, count, datatype, root, comm, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: sendbuf, recvbuf
    INTEGER, INTENT(IN) :: count, root
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: datatype
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_gather_f   ( sendbuf, recvbuf, count, datatype, root, comm, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Gather

  !int waitio_mpi_allgather ( const void  *sendbuf,  void  *recvbuf,  int count,  WAITIO_MPI_Datatype datatype,
  !                                 WAITIO_MPI_Comm comm);
  subroutine WAITIO_MPI_Allgather ( sendbuf, recvbuf, count, datatype, comm, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: sendbuf, recvbuf
    INTEGER, INTENT(IN) :: count
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: datatype
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_allgather_f   ( sendbuf, recvbuf, count, datatype, comm, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Allgather

  ! int waitio_mpi_scatter ( void *sendbuf, int sendcount, WAITIO_MPI_Datatype sendtype,
  !                               void *recvbuf, int recvcount, WAITIO_MPI_Datatype recvtype,
  !                               int root, WAITIO_MPI_Comm comm);
  subroutine WAITIO_MPI_Scatter ( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: sendbuf, recvbuf
    INTEGER, INTENT(IN) :: sendcount, recvcount, root
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: sendtype, recvtype
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_scatter_f ( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Scatter

  ! int waitio_mpi_scatterv ( const void  *sendbuf,  const int *sendcounts,  const int *displs,  
  !                           WAITIO_MPI_Datatype sendtype,  void  *recvbuf,  int recvcount,  
  !                           WAITIO_MPI_Datatype recvtype,  int root,  WAITIO_MPI_Comm comm);
  subroutine WAITIO_MPI_Scatterv ( sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: sendbuf, recvbuf, sendcounts, displs
    INTEGER, INTENT(IN) :: recvcount, root
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: sendtype, recvtype
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_scatter_f ( sendbuf, sendcounts, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Scatterv

  ! int waitio_mpi_gatherv ( const void  *sendbuf,  int sendcount,  WAITIO_MPI_Datatype sendtype,  
  !                          void  *recvbuf,  const int *recvcounts,  const int *displs,  
  !                          WAITIO_MPI_Datatype recvtype,  int root,  WAITIO_MPI_Comm comm);
  subroutine WAITIO_MPI_Gatherv ( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: sendbuf, recvbuf, recvcounts, displs
    INTEGER, INTENT(IN) :: sendcount, root
    TYPE(WAITIO_MPI_Datatype), INTENT(IN) :: sendtype, recvtype
    TYPE(WAITIO_MPI_Comm), INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_mpi_gatherv_f  ( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, ierr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Gatherv

  subroutine WAITIO_MPI_Barrier ( comm, ierr)
    implicit none
    TYPE(WAITIO_MPI_Comm), INTENT(OUT) :: comm
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_barrier_f  ( comm, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Barrier
  
  subroutine WAITIO_Wait (req, ierr)
    implicit none
    type(*), dimension(*), INTENT(IN) :: req
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_wait_f   (req, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_Wait

  subroutine WAITIO_Finalize ( ierr)
    implicit none
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_finalize_f   (1000)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_Finalize

  subroutine WAITIO_Create_pbgroup ( comm, pbid, ierr)
    implicit none
    TYPE(WAITIO_MPI_Comm), INTENT(OUT) :: comm
    INTEGER, INTENT(IN) :: pbid
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_create_pbgroup_f  ( comm, pbid, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_Create_pbgroup

  subroutine WAITIO_Group_rank ( comm, rank, ierr)
    implicit none
    TYPE(WAITIO_MPI_Comm), INTENT(OUT) :: comm
    INTEGER, INTENT(OUT) :: rank
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_group_rank_f  ( comm, rank, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_Group_rank
  
  subroutine WAITIO_Group_size ( comm, size, ierr)
    implicit none
    TYPE(WAITIO_MPI_Comm), INTENT(OUT) :: comm
    INTEGER, INTENT(OUT) :: size
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_group_size_f  ( comm, size, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_Group_size

  subroutine WAITIO_MPI_Comm_rank ( comm, rank, ierr)
    implicit none
    TYPE(WAITIO_MPI_Comm), INTENT(OUT) :: comm
    INTEGER, INTENT(OUT) :: rank
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_group_rank_f  ( comm, rank, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Comm_rank
  
  subroutine WAITIO_MPI_Comm_size ( comm, size, ierr)
    implicit none
    TYPE(WAITIO_MPI_Comm), INTENT(OUT) :: comm
    INTEGER, INTENT(OUT) :: size
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_group_size_f  ( comm, size, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_MPI_Comm_size

  subroutine WAITIO_PB_rank ( rank, ierr)
    implicit none
    INTEGER, INTENT(OUT) :: rank
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_pb_rank_f  ( rank, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_PB_rank
  
  subroutine WAITIO_PB_size ( size, ierr)
    implicit none
    INTEGER, INTENT(OUT) :: size
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr
    integer :: cerr

    call waitio_pb_size_f  ( size, cerr)
    if(present(ierr)) ierr = cerr;

  end subroutine WAITIO_PB_size

end module waitio_mpif90
