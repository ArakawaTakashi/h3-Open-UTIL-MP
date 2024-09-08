int NUM_PB  ; 
int PB_ID      ; 
int *pb_nprocs ;
int *king_rank ;


struct _comm_t {
  int             comm_id ; 
  waitio_group_t  waitio_comm ;
  MPI_Comm        mpi_comm ;
  int             is_univ ; // num_of_world == 1 or > 1
  int             my_size_universe ;
  int             my_rank_universe ;
  int             my_size_world    ;
  int             my_rank_world    ; 
  int             leader_rank ; //universal rank of my leader
  int             num_of_world ; // number of my world
  int            *world_size   ; // size(nummber of rank) of each world
  int            *world_king   ; // universal rank of each world king
  int             my_world_id  ; // 0 <= id < num_of_world
  struct _comm_t *next_ptr ; 
  struct _comm_t *parent_ptr ; 
} ; 

typedef struct _comm_t comm_t ;


comm_t *the_universe ;

comm_t *start_comm   ; // start_comm == the_universe
comm_t *last_comm    ; //
comm_t *current_comm ;

int init_comm(void) ;

int print_comm_info(comm_t *comm)  ; 

comm_t *get_comm(int comm_id)  ; 

comm_t *get_last_comm(void) ;

waitio_group_t get_waitio_comm(comm_t *comm)  ;

MPI_Comm get_mpi_comm(comm_t *comm)  ; 

int get_my_rank_universe(comm_t *comm) ;

int get_my_size_universe(comm_t *comm) ;

int get_my_rank_world(comm_t *comm) ;

int get_my_size_world(comm_t *comm) ;

int get_my_world_id(comm_t *comm) ;

int get_num_of_world(comm_t *comm) ;

int get_world_id_from_rank(comm_t *comm, int rank)  ; 

int get_world_rank_from_rank(comm_t *comm, int urank)  ; 

int get_world_size(comm_t *comm, int world_id) ;

int get_king_rank(comm_t *comm, int world_id) ;

int get_universal_rank_from_rank(comm_t *comm, int world_id, int wrank)  ; 


int split_comm(comm_t comm, int color, int key, comm_t * new_comm) ;

int allgather_int(comm_t comm, int num, int * int_array)  ; 

