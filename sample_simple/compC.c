#include "h3ouc.h"
#include <stdio.h>

int main(void) {
  int target_data ;
  int my_data ; 
  int send_data[5] ; 
  int recv_data[5] ;
  int array_size = 5 ; 
  int recv_array[10] ;
  
  my_data = 100 ;
  
  h3ouc_init("compC", "coupling.conf") ;
  h3ouc_recv_int("compA", & target_data) ;
  h3ouc_send_int("compA", my_data) ; 
  printf(" compC target_data = %d\n", target_data) ; 
  h3ouc_recv_int_array("compD", recv_data, array_size) ; 
  for (int i = 0 ; i < array_size ; i++) {
    printf("%d\n", recv_data[i]) ;
    send_data[i] = i ; 
  }

  h3ouc_send_int_array("compD", send_data, array_size) ;

  // h3ouc_irecv_model_int("compA", 0, recv_array, 10) ;
  //h3ouc_irecv_waitall() ;
  
  h3ouc_end() ;
  return 0 ; 
}
