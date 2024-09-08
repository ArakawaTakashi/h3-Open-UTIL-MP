#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "mpi.h"
#include "waitio.h"
#include "impi_comm.h"
#include "impi.h"


struct _impi_req_t {
  waitio_req_t wreq ;
  int         index ;    
  bool     use_flag ;
} ;


typedef struct _impi_req_t impi_req_t ;

const int MAX_REQ = 1000 ;

impi_req_t *req_array ;


//------------------------------------------------------------------------------------------------

int iMPI_Init(void){
  init_comm() ;

  iMPI_COMM_UNIVERSE = the_universe->waitio_comm ;
  iMPI_COMM_WORLD    = the_universe->comm_id ; 

  req_array = (impi_req_t *)malloc(sizeof(impi_req_t)*MAX_REQ) ;

  impi_req_t *req ; 
  for (int i = 0 ; i<MAX_REQ ; i++) {
    req = (req_array + i) ; 
    req->use_flag = false ;
    req->index    = i ;
  } ;
    
  return 0 ;
  
}

//------------------------------------------------------------------------------------------------

impi_req_t * get_new_request(void) {
  impi_req_t *now_req ; 
  for (int i = 0 ; i < MAX_REQ ; i++) {
    now_req = (req_array + i) ;
    if (now_req->use_flag == false) {
      now_req->use_flag = true ; 
      return now_req ;
    }
  } ;
  fprintf(stderr, "ERROR !!!!!, [impi.c] get_new_request, request size exeeded MAX_REQ\n") ; 
}

//------------------------------------------------------------------------------------------------

int iMPI_Finalize(void) {
  free(req_array) ; 
  waitio_finalize() ; 
  MPI_Finalize() ;
  return 0 ; 
}


//------------------------------------------------------------------------------------------------

int iMPI_Get_comm_world(int *world_comm) {
  * world_comm = iMPI_COMM_WORLD ;
  return iMPI_COMM_WORLD ; 
}

//------------------------------------------------------------------------------------------------

int iMPI_Comm_split(int comm_id, int color, int key, int *new_comm_id) {
  comm_t * current_comm ;

  current_comm = get_comm(comm_id) ;

  comm_t * new_comm = (comm_t *)malloc(sizeof(comm_t)) ;

  split_comm(*current_comm, color, key, new_comm) ;

  *new_comm_id = new_comm->comm_id ;
  
  return *new_comm_id ; 
}  
		    
//------------------------------------------------------------------------------------------------

int iMPI_Comm_rank(int comm_id, int *rank) {
  comm_t * current_comm ;

  current_comm = get_comm(comm_id) ;

  *rank = current_comm->my_rank_universe ; 
  return *rank ; 
}

//------------------------------------------------------------------------------------------------

int iMPI_Comm_size(int comm_id, int *size) {
  comm_t *current_comm ;

  current_comm = get_comm(comm_id) ;

  *size = current_comm->my_size_universe ; 
  return *size ; 
}

//------------------------------------------------------------------------------------------------

int iMPI_Isend(void *buffer, int count, MPI_Datatype datatype, int dest,  int tag, IMPI_Comm comm_id, IMPI_Request *req) {
  comm_t * comm = get_comm(comm_id) ;
  
  int datasize ;
  MPI_Type_size(datatype, &datasize) ;

  impi_req_t * new_req = get_new_request() ; 
  waitio_isend(get_waitio_comm(comm), dest, (char *)buffer, count*datasize, tag, &(new_req->wreq)) ;

  *req = new_req->index ; 

  return count ; 
}

//------------------------------------------------------------------------------------------------

int iMPI_Irecv(void *buffer, int count, MPI_Datatype datatype, int source, int tag, IMPI_Comm comm_id, IMPI_Request *req) {
  comm_t * comm = get_comm(comm_id) ;

  int datasize ;
  MPI_Type_size(datatype, &datasize) ;

  impi_req_t * new_req = get_new_request() ; 
  waitio_irecv(get_waitio_comm(comm), source, (char *)buffer, count*datasize, tag, &(new_req->wreq)) ;

  *req = new_req->index ;
  
  return count ; 
}

//------------------------------------------------------------------------------------------------

int iMPI_wait(IMPI_Request *req, MPI_Status *stat) {

  int current_index = *req ; 
  impi_req_t * now_req = (req_array + current_index) ; 
  waitio_wait(&(now_req->wreq)) ;
  now_req->use_flag = false ;

  return now_req->index ; 
}

//------------------------------------------------------------------------------------------------

int iMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, IMPI_Comm comm_id) {
  comm_t * comm = get_comm(comm_id) ;

  //fprintf(stderr, "iMPI_Bcast count = %d, datatype = %d, root = %d, comm_id = %d\n", count, datatype, root, comm_id) ; 

  int datasize ;
  MPI_Type_size(datatype, &datasize) ;

  //fprintf(stderr, "iMPI_Bcast data size = %d\n", datasize) ;
  
  int root_world_id    = get_world_id_from_rank(comm, root) ; 
  int root_world_rank  = get_world_rank_from_rank(comm, root) ;
  int my_world_id      = get_my_world_id(comm) ;
  int my_rank_universe = get_my_rank_universe(comm) ; 
  int my_rank_world    = get_my_rank_world(comm) ; 
  int num_of_world     = get_num_of_world(comm) ;
  //printf("num_of_world = %d\n", num_of_world) ;
  //printf("iMPI_Bcast, %d %d %d %d\n", my_rank_universe, my_rank_world, root_world_id, my_world_id) ;

  MPI_Request req ;
  MPI_Status  stat ;
  waitio_req_t wreq ;

  if (comm->is_univ) { // inter world communication
    if (root_world_id == my_world_id) { // send data to other world
      if (root_world_rank != 0) { // set data to the world king (rank_world == 0)
	if (my_rank_universe == root) {
	  MPI_Isend(buffer, count, datatype, 0, 0, get_mpi_comm(comm), &req) ;
	  MPI_Wait(&req, &stat) ; 
	}
	if (my_rank_world == 0) { // I'm a king
	  MPI_Irecv(buffer, count, datatype, root_world_rank, 0, get_mpi_comm(comm), &req) ;
	  MPI_Wait(&req, &stat) ; 
	}
      }
      for (int i = 0 ; i < num_of_world ; i++) {
        if (i != root_world_id) {
          if (my_rank_world == 0) {
            waitio_isend(get_waitio_comm(comm), get_king_rank(comm, i), (char *)buffer, count*datasize, 0, &wreq) ;
	    waitio_wait(&wreq) ;
          }
      	}
      }
    } else { // recv data from root world
      if (my_rank_world == 0) { // world king
        waitio_irecv(get_waitio_comm(comm), get_king_rank(comm, root_world_id), (char *)buffer, count*datasize, 0, &wreq) ;
        waitio_wait(&wreq) ;
      }
    }
  }

  int ret = MPI_Bcast(buffer, count, datatype, 0, get_mpi_comm(comm)) ; 

  return ret ; 
}

//------------------------------------------------------------------------------------------------

int iMPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
		int root, IMPI_Comm comm_id) {

  comm_t *comm = get_comm(comm_id) ;
  
  int datasize ;
  MPI_Type_size(datatype, &datasize) ;

  int root_world_id    = get_world_id_from_rank(comm, root) ; 
  int root_world_rank  = get_world_rank_from_rank(comm, root) ;
  int my_world_id      = get_my_world_id(comm) ;
  int my_rank_universe = get_my_rank_universe(comm) ; 
  int my_rank_world    = get_my_rank_world(comm) ; 
  int num_of_world     = get_num_of_world(comm) ;
  //printf("num_of_world = %d\n", num_of_world) ;
  //printf("iMPI_Reduce, %d %d %d %d\n", my_rank_universe, my_rank_world, root_world_id, my_world_id) ;

  MPI_Request req ;
  MPI_Status  stat ;
  waitio_req_t wreq ;
  int ret  = 0 ;
  
  if (comm->is_univ){
    char * my_world_buffer ; 

    if (my_rank_world == 0) my_world_buffer = (char *) malloc(count*datasize) ;  // intra world reduce

    ret = MPI_Reduce(sendbuf, my_world_buffer, count, datatype, op, 0, get_mpi_comm(comm)) ; // intra world reduce

    char *all_world_buffer ;
    
    if (my_rank_world == 0) { // inter world reduce
      if (my_world_id == root_world_id) { // inter world king

	all_world_buffer = (char *) malloc(count*datasize*num_of_world) ;

	for (int i = 0 ; i < num_of_world ; i++) { //recv from other world
	  if (i == my_world_id) {
	    memcpy(all_world_buffer + i*count*datasize, my_world_buffer, count*datasize) ; 
	  } else {
	    waitio_irecv(get_waitio_comm(comm), get_king_rank(comm, i), (all_world_buffer + i*count*datasize),
	    		 count*datasize, 0, &wreq) ;
	    waitio_wait(&wreq) ;
	  }
	}

	if (datatype == MPI_INT) {
	  Reduce_array_int((int *)all_world_buffer, count, num_of_world, op) ;
	} else if (datatype == MPI_LONG) {
	  Reduce_array_long((long *)all_world_buffer, count, num_of_world, op) ;
        } else if (datatype == MPI_FLOAT) {
	  Reduce_array_float((float *)all_world_buffer, count, num_of_world, op) ;
	} else if (datatype == MPI_DOUBLE) {
	  Reduce_array_double((double *)all_world_buffer, count, num_of_world, op) ;
	} else {
	  fprintf(stderr, "iMPI_Reduce, datatype invalid\n") ;
	  exit(0) ;
	}

      } else {
        waitio_isend(get_waitio_comm(comm), get_king_rank(comm, root_world_id), my_world_buffer,
		     count*datasize, 0, &wreq) ;
	waitio_wait(&wreq) ;
      }
    }

    if (my_world_id == root_world_id) {
      if (root_world_rank == 0) {
	if (my_rank_world == 0) {
	  memcpy(recvbuf, all_world_buffer, count*datasize) ; 
	}
      } else {
	if (my_rank_world == 0) {
	  MPI_Isend(all_world_buffer, count, datatype, root_world_rank, 0, get_mpi_comm(comm), &req) ;
	  MPI_Wait(&req, &stat) ; 
	} else if (my_rank_world == root_world_rank) {
	  MPI_Irecv(recvbuf, count, datatype, 0, 0, get_mpi_comm(comm), &req) ;
	  MPI_Wait(&req, &stat) ; 
	}
      }
    }

    if (my_rank_world == 0) {
      free(my_world_buffer) ;
      if (my_world_id == root_world_id) {
	free(all_world_buffer) ;
      }
    }

  } else {
    ret = MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, get_mpi_comm(comm)) ;
  }
    
  return ret ;
  
}

//------------------------------------------------------------------------------------------------

int Reduce_array_int(int * array, int nx, int ny, MPI_Op op) {

  if (op == MPI_MAX) {
    for (int i = 0 ; i < nx ; i++) {
      int max = INT_MIN ;
      for (int j = 0 ; j < ny ; j++) {
	if (*(array + i + nx*j) > max) max = *(array + i + nx*j) ;
      }
      *(array + i) = max ;
    }
  } else if (op == MPI_MIN) {
    for (int i = 0 ; i < nx ; i++) {
      int min = INT_MAX ;
      for (int j = 0 ; j < ny ; j++) {
	if (*(array + i + nx*j) < min) min = *(array + i + nx*j) ;
      }
      *(array + i) = min ;
    }
  } else if (op == MPI_SUM) {
    for (int i = 0 ; i < nx ; i++) {
      int sum = 0 ; 
      for (int j = 0 ; j < ny ; j++) {
	sum += *(array + i + nx*j) ; 
      }
      *(array + i) = sum ;
    }
  } else {
    fprintf(stderr, "Reduce_array_int, MPI_Op invalid\n") ; 
    exit(0) ;
  }

  return 0 ; 
  
}

//------------------------------------------------------------------------------------------------

int Reduce_array_long(long * array, int nx, int ny, MPI_Op op) {

  if (op == MPI_MAX) {
    for (int i = 0 ; i < nx ; i++) {
      long max = LONG_MIN ;
      for (int j = 0 ; j < ny ; j++) {
	if (*(array + i + nx*j) > max) max = *(array + i + nx*j) ;
      }
      *(array + i) = max ;
    }

  } else if (op == MPI_MIN) {
    for (int i = 0 ; i < nx ; i++) {
      long min = LONG_MAX ;
      for (int j = 0 ; j < ny ; j++) {
	if (*(array + i + nx*j) < min) min = *(array + i + nx*j) ;
      }
      *(array + i) = min ;
    }
  } else if (op == MPI_SUM) {
    for (int i = 0 ; i < nx ; i++) {
      long sum = 0 ; 
      for (int j = 0 ; j < ny ; j++) {
	sum += *(array + i + nx*j) ; 
      }
      *(array + i) = sum ;
    }
  } else {
    fprintf(stderr, "Reduce_array_long, MPI_Op invalid\n") ; 
    exit(0) ;
  }

  return 0 ; 
  
}

//------------------------------------------------------------------------------------------------

int Reduce_array_float(float * array, int nx, int ny, MPI_Op op) {
  if (op == MPI_MAX) {
    for (int i = 0 ; i < nx ; i++) {
      float max = FLT_MIN ;
      for (int j = 0 ; j < ny ; j++) {
	if (*(array + i + nx*j) > max) max = *(array + i + nx*j) ;
      }
      *(array + i) = max ;
    }
  } else if (op == MPI_MIN) {
    for (int i = 0 ; i < nx ; i++) {
      float min = FLT_MAX ;
      for (int j = 0 ; j < ny ; j++) {
	if (*(array + i + nx*j) < min) min = *(array + i + nx*j) ;
      }
      *(array + i) = min ;
    }
  } else if (op == MPI_SUM) {
    for (int i = 0 ; i < nx ; i++) {
      float sum = 0 ; 
      for (int j = 0 ; j < ny ; j++) {
	sum += *(array + i + nx*j) ; 
      }
      *(array + i) = sum ;
    }
  } else {
    fprintf(stderr, "Reduce_array_float, MPI_Op invalid\n") ; 
    exit(0) ;
  }

  return 0 ; 
  
}

//------------------------------------------------------------------------------------------------

int Reduce_array_double(double * array, int nx, int ny, MPI_Op op) {
  if (op == MPI_MAX) {
    for (int i = 0 ; i < nx ; i++) {
      double max = DBL_MIN ;
      for (int j = 0 ; j < ny ; j++) {
	if (*(array + i + nx*j) > max) max = *(array + i + nx*j) ;
      }
      *(array + i) = max ;
    }
  } else if (op == MPI_MIN) {
    for (int i = 0 ; i < nx ; i++) {
      double min = DBL_MAX ;
      for (int j = 0 ; j < ny ; j++) {
	if (*(array + i + nx*j) < min) min = *(array + i + nx*j) ;
      }
      *(array + i) = min ;
    }
  } else if (op == MPI_SUM) {
    for (int i = 0 ; i < nx ; i++) {
      double sum = 0 ; 
      for (int j = 0 ; j < ny ; j++) {
	sum += *(array + i + nx*j) ; 
      }
      *(array + i) = sum ;
    }
  } else {
    fprintf(stderr, "Reduce_array_double, MPI_Op invalid\n") ; 
    exit(0) ;
  }

  return 0 ; 
  
}

//------------------------------------------------------------------------------------------------

int iMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
		   IMPI_Comm comm_id) {

  iMPI_Reduce(sendbuf, recvbuf, count, datatype, op, 0, comm_id) ; 
  iMPI_Bcast(recvbuf, count, datatype, 0, comm_id) ;
  
}

//------------------------------------------------------------------------------------------------

int iMPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, IMPI_Comm comm_id) {
  comm_t * comm = get_comm(comm_id) ;

  int datasize ;
  MPI_Type_size(sendtype, &datasize) ;

  int root_world_id    = get_world_id_from_rank(comm, root) ; 
  int root_world_rank  = get_world_rank_from_rank(comm, root) ;
  int my_world_id      = get_my_world_id(comm) ;
  int my_rank_universe = get_my_rank_universe(comm) ; 
  int my_rank_world    = get_my_rank_world(comm) ; 
  int num_of_world     = get_num_of_world(comm) ;
  //printf("num_of_world = %d\n", num_of_world) ;
  //printf("iMPI_Bcast, %d %d %d %d\n", my_rank_universe, my_rank_world, root_world_id, my_world_id) ;

  MPI_Request req ;
  MPI_Status  stat ;
  waitio_req_t wreq ;
  int ret = 0 ;
  
  if (comm->is_univ) { // inter world communication

    int recvcount_world = 0 ;
    char *recvbuf_world ;
    
    // intra world gather
    if (my_rank_world == 0) { 
      recvcount_world = sendcount*get_my_size_world(comm) ; 
      recvbuf_world = (char *)malloc(datasize*recvcount_world) ;
    }
    ret = MPI_Gather(sendbuf, sendcount, sendtype, recvbuf_world, sendcount, recvtype, 0, get_mpi_comm(comm)) ;  // gather to the world king

    // inter world gather (send/recv)

    int total_size = 0 ;
    for (int i = 0 ; i < num_of_world ; i++) {
      total_size += get_world_size(comm, i) ;
    }

    char *all_world_buffer ; 

    if (my_rank_world == 0) {
      if (my_world_id == root_world_id) { // gather to the root world

	all_world_buffer = (char *)malloc(total_size*sendcount*datasize) ; 

	int offset = 0 ; 
	for (int i = 0 ; i < num_of_world ; i++) {
	  int world_size = get_world_size(comm, i) ; 
	  if (i == my_world_id) {
	    memcpy(all_world_buffer+offset, recvbuf_world, recvcount_world*datasize) ; 
	    offset += recvcount_world*datasize ; 
	  } else {
	    waitio_irecv(get_waitio_comm(comm), get_king_rank(comm, i), (char *)(all_world_buffer + offset),
	    		 sendcount*world_size*datasize, 0, &wreq) ;
	    waitio_wait(&wreq) ;
	    offset += sendcount*world_size*datasize ; 
	  }
	}
      } else {
	waitio_isend(get_waitio_comm(comm), get_king_rank(comm, root_world_id), recvbuf_world,
		     recvcount_world*datasize, 0, &wreq) ;
	waitio_wait(&wreq) ;
      }
    }

    //copy from 0 to root
    if (my_world_id == root_world_id) {
      if (root_world_rank == 0) {
	if (my_rank_world == 0) {
	  memcpy(recvbuf, all_world_buffer, total_size*sendcount*datasize) ; 
	}
      } else {
	if (my_rank_world == 0) {
	  MPI_Isend(all_world_buffer, total_size*sendcount, sendtype, root_world_rank, 0, get_mpi_comm(comm), &req) ;
	  MPI_Wait(&req, &stat) ; 
	} else if (my_rank_world == root_world_rank) {
	  MPI_Irecv(recvbuf, total_size*sendcount, sendtype, 0, 0, get_mpi_comm(comm), &req) ;
	  MPI_Wait(&req, &stat) ; 
	}
      }
    }
      
    if (my_rank_world == 0) {
      free(recvbuf_world) ; 
      if (my_world_id == root_world_id) {
	free(all_world_buffer) ;
      }
    }

  } else { // intra world communication only
    ret = MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, get_mpi_comm(comm)) ;
  }
    
  return ret ;

}

//------------------------------------------------------------------------------------------------

int iMPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
		   void *recvbuf, int recvcount, MPI_Datatype recvtype, IMPI_Comm comm_id) {
  iMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, 0, comm_id) ; 

  comm_t * comm = get_comm(comm_id) ;
  int total_size = 0 ;
  for (int i = 0 ; i < get_num_of_world(comm) ; i++) {
    total_size += get_world_size(comm, i) ;
  }
  
  iMPI_Bcast(recvbuf, sendcount*total_size, recvtype, 0, comm_id) ;
}

//------------------------------------------------------------------------------------------------

int get_pb_num_from_rank(int rank) {
  int min = 0 ;
  int max = 0 ; 
  for (int i = 0 ; i < NUM_PB ; i++) {
    max = min + *(pb_nprocs + i) - 1  ;    
    if ((min <= rank) && (rank <= max)) return i + 1 ;
    min = max + 1 ;
  }
}

/*
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
        }
    }
}
*/

/*
 * FORTRAN Interface
 *
 */
/*
void impi_init_ (IMPI_Fint  *ierr) {
  *ierr = iMPI_Init() ; 
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

void impi_bcast_ (void *buffer, int *count, IMPI_Datatype *datatype, int *root, IMPI_Comm *comm_id, int *ierr) {
  iMPI_Bcast(buffer, *count, impi_dtype(*datatype), *root, *comm_id) ;
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_isend_ (void *buffer, int *count, IMPI_Datatype *datatype, int *dest, int *tag, IMPI_Comm *comm_id, MPI_Request *req, int *ierr) {
  iMPI_Isend(buffer, *count, impi_dtype(*datatype), *dest, *tag, *comm_id, req) ;
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_irecv_ (void *buffer, int *count, IMPI_Datatype *datatype, int *source, int *tag, IMPI_Comm *comm_id, MPI_Request *req, int *ierr) {
  iMPI_Irecv(buffer, *count, impi_dtype(*datatype), *source, *tag, *comm_id, req) ;
  *ierr = 0 ;
}

//------------------------------------------------------------------------------------------------

void impi_iwait_ (MPI_Request *req, MPI_Status *stat, int *ierr) {
  iMPI_wait(req, stat) ; 
  *ierr = 0 ;
}

*/
