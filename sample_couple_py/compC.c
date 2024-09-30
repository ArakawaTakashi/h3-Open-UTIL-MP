#include "h3ouc.h"
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

//int MAIN__(void) {
int main(void) {
  char * my_name = "compC" ;
  const int global_array_size = 20 ;
  int time_array[6] = {2000, 1, 1, 0, 0, 0};
  int delta_t       = 3600 ; 
  int    array_size ;
  int    ierror ;
  
  // initialize  
  h3ouc_init(my_name, "coupling.conf") ;

  // get mpi parameters

  int my_rank = h3ouc_get_my_rank() ;
  int my_size = h3ouc_get_my_size() ; 

  printf("compC my_size = %d, my_rank = %d\n", my_size, my_rank) ;

  if ((global_array_size % my_size) != 0) {
    printf("Error, mod(global_array_size, my_size) must be 0\n") ;
    exit(ierror) ;
  }

  array_size = global_array_size/my_size ;
  
  int    my_grid[array_size] ; 
  int    my_array_int[array_size] ;
  float  my_array_float[array_size] ;
  double my_array_double[array_size] ;
  int    target_array_int[array_size] ;
  float  target_array_float[array_size] ;
  double target_array_double[array_size] ;
  
  for (int i = 0 ; i < array_size ; i++) {
    my_grid[i] = array_size * my_rank + i + 1 ; 
  }

  h3ouc_def_grid(my_grid, array_size, my_name, "compC_grid1", 20) ;
  h3ouc_end_grid_def() ; 

  h3ouc_set_interpolation_table_no_index(my_name, "compA", "compA_grid1", "compC", "compC_grid1", 1) ;

  h3ouc_set_interpolation_table_no_index(my_name, "compC", "compC_grid1", "compA", "compA_grid1", 1) ;

  h3ouc_init_time(time_array) ;

  // set my array
  for (int i = 0 ; i < array_size ; i++) {
    my_array_int[i] = 30000 + my_rank * 100 + i +1 ; 
    my_array_float[i] = my_array_int[i] ;
    my_array_double[i] = my_array_int[i] ;
  }

  h3ouc_put_data_1d("compC_var1", my_array_double) ;

  int is_get_ok ;
  
  for (int t = 0 ; t < 24 ; t++) {
    h3ouc_set_time(my_name, time_array, delta_t) ;

    h3ouc_get_data_1d("compA_var2", target_array_double, & is_get_ok) ;

    h3ouc_put_data_1d("compC_var1", my_array_double) ;

    h3ouc_inc_calendar(time_array, delta_t) ; 
    
  }

  h3ouc_coupling_end(time_array, 1) ;

  return 0 ; 


}


