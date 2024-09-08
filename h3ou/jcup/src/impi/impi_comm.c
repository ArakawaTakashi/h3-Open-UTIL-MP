#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "waitio.h"
#include "impi_comm.h"

waitio_group_t WAITIO_COMM_UNIVERSE ;
MPI_Comm WAITIO_COMM_WORLD    ;

//------------------------------------------------------------------------------------------------

int init_comm(void) {
  waitio_filter_func_t func[4]= {truef, truef, truef, truef};
  int  array[4] = {1, 2, 3, 4};
  int ret, size, rank, wsize, wrank;
  int num_pb, pb_id;

  waitio_init(1000) ;

  waitio_pb_size(&num_pb);
  waitio_pb_rank(&pb_id);

  //fprintf(stderr, "num_pb = %d, pb_id = %d\n", num_pb, pb_id) ;
  
  WAITIO_COMM_UNIVERSE = waitio_create_group(0, func, array);
  if(WAITIO_COMM_UNIVERSE == NULL) {
    ret = -1 ; 
    fprintf(stderr, "waitio_create_group failed code %d\n", ret);
    MPI_Finalize();
    exit(ret);
  }
  waitio_group_size(WAITIO_COMM_UNIVERSE, &wsize);
  waitio_group_rank(WAITIO_COMM_UNIVERSE, &wrank);

  NUM_PB = num_pb ;
  PB_ID  = pb_id ;

  pb_nprocs = (int *)malloc(sizeof(int)*NUM_PB) ;
  king_rank = (int *)malloc(sizeof(int)*NUM_PB) ; 
  waitio_get_nprocs(pb_nprocs) ;

  int sum = 0 ;
  for (int i = 0 ; i < NUM_PB; i++) {
    *(king_rank + i) = sum ;
    sum += *(pb_nprocs + i) ;
  } ; 

  fprintf(stderr, "sum = %d\n", sum) ; 

  //set the_universe comm
  the_universe = (comm_t *)malloc(sizeof(comm_t)) ;
  the_universe->comm_id           = 0 ; 
  the_universe->waitio_comm       = WAITIO_COMM_UNIVERSE; 
  the_universe->mpi_comm          = MPI_COMM_WORLD ; 
  the_universe->is_univ           = 1 ;
  MPI_Comm_rank(the_universe->mpi_comm, &(the_universe->my_rank_world)) ; 
  MPI_Comm_size(the_universe->mpi_comm, &(the_universe->my_size_world)) ;

  waitio_group_rank(WAITIO_COMM_UNIVERSE, &the_universe->my_rank_universe);
  waitio_group_size(WAITIO_COMM_UNIVERSE, &the_universe->my_size_universe);
#ifdef ORIG
  //the_universe->my_rank_universe  = *(king_rank + pb_id - 1) + the_universe->my_rank_world ; 
  //the_universe->my_size_universe  = sum ;
#endif
  fprintf(stderr, "wrank=%d, wsize=%d, cal-wrank=%d, cal-wsize=%d\n",
	  the_universe->my_rank_universe,  the_universe->my_size_universe,
	  *(king_rank + pb_id - 1) + the_universe->my_rank_world, sum  );

  the_universe->leader_rank = 0 ;
  the_universe->num_of_world = num_pb ;
  the_universe->world_size = (int *)malloc(sizeof(int)*num_pb) ;
  for (int i = 0 ; i < num_pb ; i++) {
    *(the_universe->world_size + i) = *(pb_nprocs + i) ;
  }
  the_universe->world_king = (int *)malloc(sizeof(int)*num_pb) ;
  int king_num = 0 ; 
  for (int i = 0 ; i < num_pb ; i++) {
    *(the_universe->world_king + i) = king_num ;
    king_num += *(the_universe->world_size + i) ;
  }
  the_universe->my_world_id = get_world_id_from_rank(the_universe, the_universe->my_rank_universe) ; 
  the_universe->parent_ptr = NULL ;
  the_universe->next_ptr   = NULL ; 

  start_comm   = the_universe ; 
  last_comm    = the_universe ;
  current_comm = the_universe ; 

  //print_comm_info(start_comm) ;

  return 0 ; 
}
  
//------------------------------------------------------------------------------------------------

int print_comm_info(comm_t *comm) {
  fprintf(stderr, "--------  comm info  --------\n") ; 
  fprintf(stderr, "  comm_id     = %d\n", comm->comm_id) ; 
  fprintf(stderr, "  watio_comm  = %d\n", comm->waitio_comm) ; 
  fprintf(stderr, "  mpi_comm    = %d\n", comm->mpi_comm) ; 
  fprintf(stderr, "  is_univ     = %d\n", comm->is_univ) ; 
  fprintf(stderr, "  size_univ   = %d\n", comm->my_size_universe) ; 
  fprintf(stderr, "  rank_univ   = %d\n", comm->my_rank_universe) ; 
  fprintf(stderr, "  my  world   = %d\n", comm->my_world_id) ; 
  fprintf(stderr, "  size_world  = %d\n", comm->my_size_world) ; 
  fprintf(stderr, "  rank_world  = %d\n", comm->my_rank_world) ; 
  fprintf(stderr, "  num world   = %d\n", comm->num_of_world) ;
  fprintf(stderr, "  world size  = ") ; 
  for (int i = 0 ; i < comm->num_of_world ; i++) {
    fprintf(stderr, "%d ", *(comm->world_size + i)) ; 
  }
  fprintf(stderr, "\n") ; 
  fprintf(stderr, "  world king  = ") ; 
  for (int i = 0 ; i < comm->num_of_world ; i++) {
    fprintf(stderr, "%d ", *(comm->world_king + i)) ; 
  }
  fprintf(stderr, "\n") ;

  return 0 ;
  
}

//------------------------------------------------------------------------------------------------

int split_comm(comm_t comm, int color, int key, comm_t * new_comm) {

  new_comm->parent_ptr = &comm ;
  
  int * color_array = (int *)malloc(sizeof(int)*comm.my_size_universe) ; 

  allgather_int(comm, color, color_array) ;

  new_comm->my_size_universe = 0 ;
  new_comm->my_rank_universe = 0 ;

  int num ; 

  fprintf(stderr, "color =  %d\n", color) ; 
  for (int i = 0 ; i < comm.my_size_universe ; i++) {
    num = *(color_array + i) ; 

    fprintf(stderr, "color_array = %d\n", num) ;
    
    if (num == color) {
      if (i == comm.my_rank_universe) {
    	new_comm->my_rank_universe = new_comm->my_size_universe ;
      }
      ++new_comm->my_size_universe ; 
    }
  }
  
  int counter = 0 ; 
  int * world_flag = (int *)malloc(sizeof(int)*comm.num_of_world) ;
  for (int j = 0 ; j < comm.num_of_world ; j++) {
    *(world_flag + j) = 0 ; 
    for (int i = 0 ; i < *(comm.world_size + j) ; i++) {
      num = *(color_array + counter) ;
      if (num == color) *(world_flag + j) += 1 ; 
      counter++ ; 
    }
  }
  
  counter = 0 ;
  for (int j = 0 ; j < comm.num_of_world ; j++) {
    if (*(world_flag + j) > 0) counter++ ;
  }

  new_comm->num_of_world = counter ; 
  new_comm->world_size = (int *)malloc(sizeof(int) * counter) ;
 
  counter = 0 ;
  for (int j = 0 ; j < comm.num_of_world ; j++) {
    if (*(world_flag + j) > 0) {
      *(new_comm->world_size + counter) = *(world_flag + j) ; 
      counter++ ;
    }
  }

  new_comm->world_king = (int *)malloc(sizeof(int) * new_comm->num_of_world) ;
  
  counter = 0 ;
  for (int i = 0 ; i < new_comm->num_of_world ; i++) {
    *(new_comm->world_king + i) = counter ;
    counter += *(new_comm->world_size + i) ;
  }
  
  free(world_flag) ;
  free(color_array) ;
  
  new_comm->is_univ = (new_comm->num_of_world > 1) ;
  
  new_comm->my_world_id = get_world_id_from_rank(new_comm, new_comm->my_rank_universe) ;
  
  MPI_Comm_split(comm.mpi_comm, color, key, &(new_comm->mpi_comm)) ;
  MPI_Comm_rank(new_comm->mpi_comm, &(new_comm->my_rank_world)) ;
  MPI_Comm_size(new_comm->mpi_comm, &(new_comm->my_size_world)) ;
  
  last_comm = get_last_comm() ;
  new_comm->comm_id = last_comm->comm_id + 1 ;
  last_comm->next_ptr = new_comm ;
  
  print_comm_info(new_comm) ;

  return 0 ; 
}

//------------------------------------------------------------------------------------------------

//int add_comm(comm_t * start_comm, int new_comm, int is_univ) {
//  comm_t *current_ptr ;
//  current_ptr = start_comm ;
//  while (current_ptr) {
//    current_ptr = current_ptr->next_ptr ;
//  }
  
//  return 0 ;
  
//}

//------------------------------------------------------------------------------------------------

comm_t *get_comm(int comm_id) {
  comm_t * current_ptr ;
  current_ptr = the_universe ; 
  
  while (current_ptr != NULL) {
    if (current_ptr->comm_id == comm_id) return current_ptr ; 
    current_ptr = current_ptr->next_ptr ; 
  } ; 

  return NULL; 
}

//------------------------------------------------------------------------------------------------

comm_t *get_last_comm(void) {
  comm_t * current_ptr ;
  current_ptr = the_universe ; 
  
  while (current_ptr != NULL) {
    if (current_ptr->next_ptr == NULL) return current_ptr ; 
    current_ptr = current_ptr->next_ptr ; 
  } ; 

  return NULL ; 
}

//------------------------------------------------------------------------------------------------

waitio_group_t get_waitio_comm(comm_t *comm) {
  return comm->waitio_comm ;
}

//------------------------------------------------------------------------------------------------

MPI_Comm get_mpi_comm(comm_t *comm) {
  return comm->mpi_comm ;
}

//------------------------------------------------------------------------------------------------

int get_my_rank_universe(comm_t *comm) {
  return comm->my_rank_universe ;
}

//------------------------------------------------------------------------------------------------

int get_my_size_universe(comm_t *comm) {
  return comm->my_size_universe ;
}

//------------------------------------------------------------------------------------------------

int get_my_rank_world(comm_t *comm) {
  return comm->my_rank_world ;
}

//------------------------------------------------------------------------------------------------

int get_my_size_world(comm_t *comm) {
  return comm->my_size_world ;
}

//------------------------------------------------------------------------------------------------

int get_world_id_from_rank(comm_t *comm, int urank) {
  int min = 0 ;
  int max = 0 ; 
  for (int i = 0 ; i < comm->num_of_world ; i++) {
    max = min + *(comm->world_size + i) - 1  ;    
    if ((min <= urank) && (urank <= max)) return i ;
    min = max + 1 ;
  }

  return 0 ; 
}

//------------------------------------------------------------------------------------------------

int get_my_world_id(comm_t *comm) {
  return comm->my_world_id ;
}

//------------------------------------------------------------------------------------------------

int get_num_of_world(comm_t *comm) {
  return comm->num_of_world ;
}

//------------------------------------------------------------------------------------------------

int get_world_rank_from_rank(comm_t *comm, int urank) {
  return urank - *(comm->world_king + get_world_id_from_rank(comm, urank)) ; 
}

//------------------------------------------------------------------------------------------------

int get_world_size(comm_t *comm, int world_id) {
  return *(comm->world_size + world_id) ;
}

//------------------------------------------------------------------------------------------------

int get_king_rank(comm_t *comm, int world_id)  {
  return *(comm->world_king + world_id) ;
}

//------------------------------------------------------------------------------------------------

int get_universal_rank_from_rank(comm_t *comm, int world_id, int wrank) {
  return *(comm->world_king + world_id) + wrank ; 
}

//------------------------------------------------------------------------------------------------

int allgather_int(comm_t comm, int num, int * int_array) {
    int *world_array ; 

    if (comm.my_rank_world == 0) world_array = (int *)malloc(sizeof(int)*comm.my_size_world) ;
    
    MPI_Gather(&num, 1, MPI_INT, world_array, 1, MPI_INT, 0, comm.mpi_comm) ; // gather num to world root

    if (comm.my_rank_world == 0) { // world leader communication

      waitio_req_t req ; 

      if (comm.my_rank_universe == 0) { // recv world leader array
	for (int i = 0 ; i < *(comm.world_size + 0) ; i++) { // set my world array
	  *(int_array + i) = *(world_array + i) ;
	}
        int pos = 0 ; 
        for (int i = 1 ; i < comm.num_of_world ; i++) {
	  pos += *(comm.world_size + i - 1) ; 
	  int array_size = *(comm.world_size + i) ; 
	  waitio_irecv(comm.waitio_comm, pos, (char *)(int_array + pos), array_size * sizeof(int), 0, &req) ;  
	  waitio_wait(&req) ; 
        }
      } else { // send world leader array to universal leader
        int array_size = *(comm.world_size + comm.my_world_id) ; 
        waitio_isend(comm.waitio_comm, 0, (char *)world_array, array_size*sizeof(int), 0, &req) ;
        waitio_wait(&req) ; 
      }

      if (comm.my_rank_universe == 0) { // send gathered array to world leaders
        int pos = 0 ; 
	for (int i = 1 ; i < comm.num_of_world ; i++) {
	  pos += *(comm.world_size + i - 1) ; 
	  waitio_isend(comm.waitio_comm, pos, (char *)int_array, comm.my_size_universe * sizeof(int), 0, &req) ;
	  waitio_wait(&req) ;
	}
      } else { // recv array from universal leader
	waitio_irecv(comm.waitio_comm, 0, (char *)int_array, comm.my_size_universe * sizeof(int),0, &req) ;  
	waitio_wait(&req) ; 
      }

      free(world_array) ; 
    }

    MPI_Bcast(int_array, comm.my_size_universe, MPI_INT, 0, comm.mpi_comm) ;

    return 0 ; 
}



//------------------------------------------------------------------------------------------------

int bcast_array(void *buffer, int count, MPI_Datatype datatype, int root, comm_t comm) {
  if (comm.is_univ) {
  } else {
    MPI_Bcast(buffer, count, datatype, root, comm.mpi_comm) ; 
  }

  return 0 ;
  
}

//------------------------------------------------------------------------------------------------




