#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "mpi.h"
#include "waitio.h"
#include "impi_comm.h"
#include "impi.h"

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------

static inline MPI_Datatype impi_dtype(IMPI_Datatype datatype) {
    if(datatype) {
      switch (datatype) {
        case IMPI_CHAR:
	  return MPI_CHAR ; 
        case IMPI_INT:
	  return MPI_INT ;
	case IMPI_LONG :
	  return MPI_LONG ;
	case IMPI_FLOAT :
	  return MPI_REAL ; 
        case IMPI_DOUBLE :
	  return MPI_DOUBLE ;
        default:
	  fprintf(stderr, "ERROR !!!!, impi_dtype, MPI_Datatype error, defined data type = %d\n", datatype) ;
	  MPI_Finalize() ;
	  exit(1) ; 
      }
    }
}

//------------------------------------------------------------------------------------------------

static inline MPI_Op impi_optype(IMPI_Op op) {
    if(op) {
      switch (op) {
        case IMPI_SUM:
	  return MPI_SUM ; 
        case IMPI_MAX:
	  return MPI_MAX ;
	case IMPI_MIN :
	  return MPI_MIN ;
        default:
	  fprintf(stderr, "ERROR !!!!, impi_optype, MPI_Op error, defined op type = %d\n", op) ;
	  MPI_Finalize() ;
	  exit(1) ; 
      }
    }
}

/*
 * FORTRAN Interface
 *
 */
void impi_init_ (IMPI_Fint  *ierr) {
  *ierr = iMPI_Init() ; 
  return;
}

//------------------------------------------------------------------------------------------------

void impi_finalize_ (IMPI_Fint  *ierr) {
  *ierr = iMPI_Finalize() ; 
  return;
}

//------------------------------------------------------------------------------------------------

void impi_get_comm_world_ (int *world_comm, int *ierr) {
  iMPI_Get_comm_world(world_comm) ; 
  *ierr = 0 ; 
}

//------------------------------------------------------------------------------------------------

void impi_comm_split_ (int *comm_id, int *color, int * key, int *new_comm_id, int *ierr) {
  
  iMPI_Comm_split(*comm_id, *color, *key, new_comm_id) ;
  *ierr = 0 ;

}  
		    
//------------------------------------------------------------------------------------------------

void impi_comm_rank_ (int *comm_id, int *rank, int * ierr) {

  iMPI_Comm_rank(*comm_id, rank)  ;
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_comm_size_ (int *comm_id, int *size, int *ierr) {

  iMPI_Comm_size(*comm_id, size) ; 
  *ierr = 0 ; 

}

//------------------------------------------------------------------------------------------------

void impi_isend_ (void *buffer, int *count, IMPI_Datatype *datatype, int *dest, int *tag, IMPI_Comm *comm_id, IMPI_Request *req, int *ierr) {
  iMPI_Isend(buffer, *count, impi_dtype(*datatype), *dest, *tag, *comm_id, req) ;
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_irecv_ (void *buffer, int *count, IMPI_Datatype *datatype, int *source, int *tag, IMPI_Comm *comm_id, IMPI_Request *req, int *ierr) {
  iMPI_Irecv(buffer, *count, impi_dtype(*datatype), *source, *tag, *comm_id, req) ;
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_wait_ (IMPI_Request *req, MPI_Status *stat, int *ierr) {
  iMPI_wait(req, stat) ; 
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_bcast_ (void *buffer, int *count, IMPI_Datatype *datatype, int *root, IMPI_Comm *comm_id, int *ierr) {
  iMPI_Bcast(buffer, *count, impi_dtype(*datatype), *root, *comm_id) ;
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_reduce_ (const void *sendbuf, void *recvbuf, int *count, IMPI_Datatype *datatype, IMPI_Op *op,
		   int *root, IMPI_Comm *comm_id, int *ierr) {
  iMPI_Reduce(sendbuf, recvbuf, *count, impi_dtype(*datatype), impi_optype(*op), *root, *comm_id) ;
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_allreduce_ (const void *sendbuf, void *recvbuf, int *count, IMPI_Datatype *datatype, IMPI_Op *op,
		      IMPI_Comm *comm_id, int *ierr) {
  iMPI_Allreduce(sendbuf, recvbuf, *count, impi_dtype(*datatype), impi_optype(*op), *comm_id)  ; 
  *ierr = 0 ; 
}

//------------------------------------------------------------------------------------------------

