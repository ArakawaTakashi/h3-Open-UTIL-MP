#include "h3ouc.h"
#include <stdio.h>
#include "mpi.h"

int main(void) {
  char * my_name = "compC" ;
  int    array_size = 10 ;
  int    my_array_int[array_size] ;
  float  my_array_float[array_size] ;
  double my_array_double[array_size] ;
  int    target_array_int[array_size] ;
  float  target_array_float[array_size] ;
  double target_array_double[array_size] ;
  
  // initialize  
  h3ouc_init(my_name, "coupling.conf") ;

  // get mpi parameters

  int my_rank = h3ouc_get_my_rank() ;
  int my_size = h3ouc_get_my_size() ; 
  
  printf("compC my_size = %d, my_rank = %d\n", my_size, my_rank) ;

  // set my array
  for (int i = 0 ; i < array_size ; i++) {
    my_array_int[i] = 30000 + my_rank * 100 + i +1 ; 
    my_array_float[i] = my_array_int[i] ;
    my_array_double[i] = my_array_int[i] ;
  }

  // bcast global

  h3ouc_bcast_global_int("compA", target_array_int, array_size) ;
  h3ouc_bcast_global_float("compB", target_array_float, array_size) ;

  if (my_rank == 0) {
    printf("compC bcast data \n") ; 
    h3ouc_bcast_global_double("compC", my_array_double, array_size) ;
  } else {
    h3ouc_bcast_global_double("compC", target_array_double, array_size) ;
    for (int i = 0 ; i < array_size ; i++) {
      printf("target array = %d\n", target_array_double[i]) ;
    }
  }

  // bcast model

  // h3ouc_bcast_model_float("compC", "compD", my_array_float, array_size) ;
  //h3ouc_bcast_model_double("compD", "compC", target_array_double, array_size) ;
  //for (int i = 0 ; i < array_size ; i++) {
  //   printf("bcast model, target array = %lf\n", target_array_double[i]) ;
  //}
  
  // send/recv model

  //if (my_rank == 0) {
  //  h3ouc_send_model_float("compD", 0, my_array_float, array_size) ;
  //  h3ouc_recv_model_double("compD", 0, target_array_double, array_size) ;
  //  for (int i = 0 ; i < array_size ; i++) {
  //    printf("compC send/recv model array = %lf\n", target_array_double[i]) ;
  //  }
  //}

  // bcast local
  
  MPI_Comm my_comm = h3ouc_get_my_comm() ;

  if (my_rank == 0) {
    MPI_Bcast(my_array_int, array_size, MPI_INT, 0, my_comm) ;
  } else {
    MPI_Bcast(target_array_int, array_size, MPI_INT, 0, my_comm) ;
    for (int i = 0 ; i < array_size ; i++) {
      printf("compC bcast local array = %d\n", target_array_int[i]) ;
    }
  }

  // send/recv local

  if (my_size > 1) {
    if (my_rank == 0) {
      MPI_Send(my_array_float, array_size, MPI_FLOAT, 1, 0, my_comm) ;
    }
    if (my_rank == 1){
      MPI_Status status ; 
      MPI_Recv(target_array_float, array_size, MPI_FLOAT, 0, 0, my_comm, & status) ; 
      for (int i = 0 ; i < array_size ; i++) {
	printf("compC send/recv local array = %f\n", target_array_float[i]) ;
      }
    }
  }
  
  h3ouc_end() ;

  return 0 ; 
}
