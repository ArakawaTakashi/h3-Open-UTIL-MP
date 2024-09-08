#include "mpi.h"

void h3ouc_init(char * my_name, char * config_file_name) ; 

void h3ouc_get_mpi_parameter(char * my_name, int * my_comm, int * my_group, int * my_size, int * my_rank);

void h3ouc_def_grid(int * grid_index, int ngrid, char * comp_name, char * grid_name, int nz);

void h3ouc_end_grid_def(void) ;

void h3ouc_set_interpolation_table(char * my_name, char * send_comp_name, char * send_grid_name,
                                   char * recv_comp_name, char * recv_grid_name, int map_tag,
				   int * send_index, int * recv_index, double * coef, int ngrid) ;

void h3ouc_set_interpolation_table_no_index(char * my_name, char * send_comp_name, char * send_grid_name,
                                            char * recv_comp_name, char * recv_grid_name, int map_tag) ;

void h3ouc_init_time(int * time_array) ;

void h3ouc_set_time(char * my_name, int * time_array, int delta_t) ;

void h3ouc_put_data_1d(char * data_name, double * data) ; 

void h3ouc_put_data_25d(char * data_name, double * data) ; 

void h3ouc_get_data_1d(char * data_name, double * data, int * is_get_ok) ;

void h3ouc_get_data_25d(char * data_name, double * data, int * is_get_ok) ;

int  h3ouc_get_num_of_put_data(void) ;

void h3ouc_get_put_data_name(int data_num, char * data_name) ;

int  h3ouc_get_num_of_get_data(void) ;

void h3ouc_get_get_data_name(int data_num, char * data_name) ;

void h3ouc_send_array_int(char * my_comp_name, char * recv_comp_name, int * array, int array_size) ;

void h3ouc_recv_array_int(char * my_comp_name, char * send_comp_name, int * array, int array_size) ;

void h3ouc_send_array_double(char * my_comp_name, char * recv_comp_name, double * array, int array_size) ;

void h3ouc_recv_array_double(char * my_comp_name, char * send_comp_name, double * array, int array_size) ;

void h3ouc_coupling_end(int * time_array) ; 

void h3ouc_inc_calendar(int * time_array, int delta_t) ; 

void h3ouc_init_simple(char * comp_name) ;

void h3ouc_end(void) ; 

void h3ouc_send_int(char * target_name, int data) ;

void h3ouc_send_float(char * target_name, float data) ;

void h3ouc_send_double(char * target_name, double data) ;

void h3ouc_recv_int(char * source_name, int * data) ; 

void h3ouc_recv_float(char * source_name, float * data) ; 

void h3ouc_recv_double(char * source_name, double * data) ; 

void h3ouc_send_int_array(char * target_name, int data[], int array_size) ;

void h3ouc_send_float_array(char * target_name, float data[], int array_size) ;

void h3ouc_send_double_array(char * target_name, double data[], int array_size) ;

void h3ouc_recv_int_array(char * source_name, int * data, int array_size) ; 

void h3ouc_recv_float_array(char * source_name, float * data, int array_sizs) ; 

void h3ouc_recv_double_array(char * source_name, double * data, int array_size) ; 
