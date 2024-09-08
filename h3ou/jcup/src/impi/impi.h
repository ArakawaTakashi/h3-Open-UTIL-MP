#include <stdbool.h>
#include "mpi.h"
#include "waitio.h"

typedef int IMPI_Fint ; 
typedef int IMPI_Comm ;
typedef int IMPI_Datatype ;
typedef int IMPI_Request ;
typedef int IMPI_Op ;

#define IMPI_DEFTAG         0x00080000

#define IMPI_INT    (0x0001|IMPI_DEFTAG)
#define IMPI_LONG   (0x0002|IMPI_DEFTAG)
#define IMPI_CHAR   (0x0003|IMPI_DEFTAG)
#define IMPI_FLOAT  (0x0004|IMPI_DEFTAG)
#define IMPI_DOUBLE (0x0005|IMPI_DEFTAG)
#define IMPI_TMASK  (0xffff)

#define IMPI_SUM    (0x0001|IMPI_DEFTAG)
#define IMPI_MAX    (0x0002|IMPI_DEFTAG)
#define IMPI_MIN    (0x0003|IMPI_DEFTAG)
#define IMPI_OMASK  (0xffff)

waitio_group_t iMPI_COMM_UNIVERSE ; 
int      iMPI_COMM_WORLD ; 

int iMPI_Init(void) ;

int iMPI_Finalize(void) ;

int iMPI_Get_comm_world(int *world_comm) ; 

int iMPI_Comm_rank(int comm_id, int *rank)  ; 

int iMPI_Comm_size(int comm_id, int *size)  ; 

int iMPI_Comm_split(int comm_id, int color, int key, int *new_comm_id) ;

int iMPI_Comm_rank(int comm_id, int *rank)  ; 

int iMPI_Comm_size(int comm_id, int *size)  ; 

int iMPI_Isend(void *buffer, int count, MPI_Datatype datatype, int dest,  int tag, IMPI_Comm comm_id, IMPI_Request *req) ; 

int iMPI_Irecv(void *buffer, int count, MPI_Datatype datatype, int source, int tag, IMPI_Comm comm_id, IMPI_Request *req) ; 

int iMPI_wait(IMPI_Request *req, MPI_Status *stat)  ; 

int iMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, IMPI_Comm comm_id) ; 

int iMPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, IMPI_Comm comm_od) ;

int iMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, IMPI_Comm comm_od) ;

int iMPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, IMPI_Comm comm_id) ; 

int iMPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
		   void *recvbuf, int recvcount, MPI_Datatype recvtype, IMPI_Comm comm_id) ; 


int Reduce_array_int(int *array, int nx, int ny, MPI_Op op) ;
int Reduce_array_long(long *array, int nx, int ny, MPI_Op op) ;
int Reduce_array_float(float *array, int nx, int ny, MPI_Op op) ;
int Reduce_array_double(double *array, int nx, int ny, MPI_Op op) ;

int get_pb_num_from_rank(int rank) ; 
