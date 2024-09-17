#include "h3ouc.h"
#include <stdio.h>
#include "mpi.h"

int main(void) {
  char * my_name = "compC" ;
  int    my_comm, my_group, my_size, my_rank ;  
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
  h3ouc_get_mpi_parameter(my_name, & my_comm, & my_group, & my_size, & my_rank) ; 
  
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

  h3ouc_bcast_model_float("compC", "compD", my_array_float, array_size) ;
  h3ouc_bcast_model_double("compD", "compC", target_array_double, array_size) ;
  for (int i = 0 ; i < array_size ; i++) {
     printf("bcast_model, target array = %lf\n", target_array_double[i]) ;
  }
  
  // send/recv model

  if (my_rank == 0) {
    h3ouc_send_model_float("compD", 0, my_array_float, array_size) ;
    h3ouc_recv_model_double("compD", 0, target_array_double, array_size) ;
    for (int i = 0 ; i < array_size ; i++) {
      printf("compC send/recv model array = %lf\n", target_array_double[i]) ;
    }
  }
  
  h3ouc_end() ;

  return 0 ; 
}
