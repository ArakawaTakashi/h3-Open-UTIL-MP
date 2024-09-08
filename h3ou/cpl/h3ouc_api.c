#include "h3ouc.h"
#include "h3oup.h"
#include <string.h>

int my_grid_size  = 0 ;

//------------------------------------------------------------------------------------------------

void h3ouc_init(char *my_name, char * config_file_name) {
  int name_len = strlen(my_name) ;
  int config_len = strlen(config_file_name) ; 
  h3oup_init(my_name, &name_len, config_file_name, &config_len) ; 
}


//------------------------------------------------------------------------------------------------

void h3ouc_get_mpi_parameter(char * my_name, int * my_comm, int * my_group, int * my_size, int * my_rank) {
  int name_len = strlen(my_name) ; 
  h3oup_get_mpi_parameter(my_name, &name_len, my_comm, my_group, my_size, my_rank) ;
}


//------------------------------------------------------------------------------------------------

void h3ouc_def_grid(int * grid_index, int index_len, char * my_name, char * grid_name, int nz) {
  int comp_name_len = strlen(my_name) ;
  int grid_name_len = strlen(grid_name) ; 

  my_grid_size = index_len ;

  h3oup_def_grid(grid_index, & index_len, my_name, & comp_name_len, grid_name, & grid_name_len, & nz) ;
}

//------------------------------------------------------------------------------------------------

void h3ouc_end_grid_def(void) {
  h3oup_end_grid_def() ; 
}

//------------------------------------------------------------------------------------------------

void h3ouc_set_interpolation_table(char * my_name, char * send_comp_name, char * send_grid_name,
				   char * recv_comp_name, char * recv_grid_name, int map_tag,
				   int * send_index, int * recv_index, double * coef, int nindex) {
  int my_name_len = strlen(my_name) ;
  int send_comp_name_len = strlen(send_comp_name) ;
  int send_grid_name_len = strlen(send_grid_name) ;
  int recv_comp_name_len = strlen(recv_comp_name) ;
  int recv_grid_name_len = strlen(recv_grid_name) ;

  h3oup_set_interpolation_table(my_name, & my_name_len,
			       send_comp_name, & send_comp_name_len,
			       send_grid_name, & send_grid_name_len,
			       recv_comp_name, & recv_comp_name_len,
			       recv_grid_name, & recv_grid_name_len,
				& map_tag,
				send_index, recv_index, coef, & nindex) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_set_interpolation_table_no_index(char * my_name, char * send_comp_name, char * send_grid_name,
				            char * recv_comp_name, char * recv_grid_name, int map_tag) {
  int my_name_len = strlen(my_name) ;
  int send_comp_name_len = strlen(send_comp_name) ;
  int send_grid_name_len = strlen(send_grid_name) ;
  int recv_comp_name_len = strlen(recv_comp_name) ;
  int recv_grid_name_len = strlen(recv_grid_name) ;

  h3oup_set_interpolation_table_no_index(my_name, & my_name_len,
					 send_comp_name, & send_comp_name_len,
					 send_grid_name, & send_grid_name_len,
					 recv_comp_name, & recv_comp_name_len,
					 recv_grid_name, & recv_grid_name_len,
					 & map_tag) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_init_time(int * time_array) {
  int array_size = 6 ;

  h3oup_init_time(time_array, & array_size) ; 
}

//------------------------------------------------------------------------------------------------

void h3ouc_set_time(char * my_name, int * time_array, int delta_t) {
  int my_name_len = strlen(my_name) ;
  int time_array_size = 6 ;
  
  h3oup_set_time(my_name, & my_name_len, time_array, & time_array_size, & delta_t) ;
}

//------------------------------------------------------------------------------------------------

void h3ouc_put_data_1d(char * data_name, double * data) {
  int name_len = strlen(data_name) ; 
    
  h3oup_put_data_1d(data_name, & name_len, data, & my_grid_size) ;
}

//------------------------------------------------------------------------------------------------

void h3ouc_put_data_25d(char * data_name, double * data){
  int name_len = strlen(data_name) ; 
  int vlayer ; 

  h3oup_get_vlayer(data_name, & name_len, & vlayer) ;
  
  h3oup_put_data_25d(data_name, & name_len, data, & my_grid_size, & vlayer) ;
}

//------------------------------------------------------------------------------------------------

void h3ouc_get_data_1d(char * data_name, double * data, int * is_get_ok) {
  int name_len = strlen(data_name) ;

  h3oup_get_data_1d(data_name, & name_len, data, & my_grid_size, is_get_ok) ;
}

//------------------------------------------------------------------------------------------------

void h3ouc_get_data_25d(char * data_name, double * data, int * is_get_ok) {
  int name_len = strlen(data_name) ;
  int vlayer ;

  h3oup_get_vlayer(data_name, & name_len, & vlayer) ;
  
  h3oup_get_data_25d(data_name, & name_len, data, & my_grid_size, & vlayer, is_get_ok) ;
}

//------------------------------------------------------------------------------------------------

int h3ouc_get_num_of_put_data(void) {
  int num_of_data ;

  h3oup_get_num_of_put_data(& num_of_data) ;

  return num_of_data ; 

}

//------------------------------------------------------------------------------------------------

void h3ouc_get_put_data_name(int data_num, char * data_name) {
  int  dnum ;
  int  name_len = strlen(data_name) ; 

  dnum = data_num + 1 ;
  
  h3oup_get_put_data_name(& dnum, data_name, & name_len) ; 

}

//------------------------------------------------------------------------------------------------

int h3ouc_get_num_of_get_data(void) {
  int num_of_data ;

  h3oup_get_num_of_get_data(& num_of_data) ;

  return num_of_data ; 

}

//------------------------------------------------------------------------------------------------

void h3ouc_get_get_data_name(int data_num, char * data_name) {
  int  dnum ;
  int  name_len = strlen(data_name) ; 

  dnum = data_num + 1 ; 

  h3oup_get_get_data_name(& dnum, data_name, & name_len) ; 

}

//------------------------------------------------------------------------------------------------

void h3ouc_send_array_int(char * my_comp_name, char * recv_comp_name, int * array, int array_size) {
  int my_name_len = strlen(my_comp_name) ; 
  int recv_name_len = strlen(recv_comp_name) ;

  h3oup_send_array_int(my_comp_name, & my_name_len, recv_comp_name, & recv_name_len, array, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_recv_array_int(char * my_comp_name, char * send_comp_name, int * array, int array_size) {
  int my_name_len = strlen(my_comp_name) ; 
  int send_name_len = strlen(send_comp_name) ;

  h3oup_recv_array_int(my_comp_name, & my_name_len, send_comp_name, & send_name_len, array, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_send_array_double(char * my_comp_name, char * recv_comp_name, double * array, int array_size) {
  int my_name_len = strlen(my_comp_name) ; 
  int recv_name_len = strlen(recv_comp_name) ;

  h3oup_send_array_double(my_comp_name, & my_name_len, recv_comp_name, & recv_name_len, array, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_recv_array_double(char * my_comp_name, char * send_comp_name, double * array, int array_size) {
  int my_name_len = strlen(my_comp_name) ; 
  int send_name_len = strlen(send_comp_name) ;

  h3oup_recv_array_double(my_comp_name, & my_name_len, send_comp_name, & send_name_len, array, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_coupling_end(int * time_array) {
  int time_array_size = 6 ;

  h3oup_coupling_end(time_array, & time_array_size) ;
}

//------------------------------------------------------------------------------------------------

void h3ouc_inc_calendar(int * time_array, int delta_t) {
  int time_array_size = 6 ;

  h3oup_inc_calendar(time_array, & time_array_size, & delta_t) ;

}

//------------------------------------------------------------------------------------------------

void h3ouc_init_simple(char * comp_name) {
  int comp_name_len ;
  char * log_level = "SILENT" ; 
  int log_level_len ;
  int stop_step ;
  int debug_mode ;

  comp_name_len = strlen(comp_name) ;
  log_level_len = strlen(log_level) ;
  stop_step = 0 ;
  debug_mode = 0 ;

  h3oup_init_simple(comp_name, & comp_name_len, log_level, & log_level_len, & stop_step, & debug_mode) ; 
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_end(void) {
  h3oup_end() ;
}

//------------------------------------------------------------------------------------------------

void h3ouc_send_int(char *target_name, int data) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_send_int_scalar(target_name, & name_len, & data) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_send_float(char *target_name, float data) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_send_real_scalar(target_name, & name_len, & data) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_send_double(char *target_name, double data) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_send_double_scalar(target_name, & name_len, & data) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_recv_int(char *target_name, int * data) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_recv_int_scalar(target_name, & name_len, data) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_recv_float(char *target_name, float * data) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_recv_real_scalar(target_name, & name_len, data) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_recv_double(char *target_name, double * data) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_recv_double_scalar(target_name, & name_len, data) ;
  
}
//------------------------------------------------------------------------------------------------

void h3ouc_send_int_array(char *target_name, int data[], int array_size) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_send_int_array(target_name, & name_len, data, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_send_float_array(char *target_name, float data[], int array_size) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_send_real_array(target_name, & name_len, data, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_send_double_array(char *target_name, double data[], int array_size) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_send_double_array(target_name, & name_len, data, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_recv_int_array(char *target_name, int * data, int array_size) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_recv_int_array(target_name, & name_len, data, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_recv_float_array(char *target_name, float * data, int array_size) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_recv_real_array(target_name, & name_len, data, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

void h3ouc_recv_double_array(char *target_name, double * data, int array_size) {
  int name_len ;

  name_len = strlen(target_name) ;

  h3oup_recv_double_array(target_name, & name_len, data, & array_size) ;
  
}

//------------------------------------------------------------------------------------------------

