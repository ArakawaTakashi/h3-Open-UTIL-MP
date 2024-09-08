void h3oup_init(char * my_name, int * name_len, char * config_file_name, int * config_len) ; 

void h3oup_get_mpi_parameter(char * my_name, int * name_len, int * my_comm, int * my_group, int * my_size, int * my_rank);

void h3oup_def_grid(int * grid_index, int * ngrid, char * comp_name, int * comp_name_len, char * grid_name, int * grid_name_len, int * nz);

void h3oup_end_grid_def(void) ;

void h3oup_set_interpolation_table(char * my_name, int * my_name_len,
				   char * send_comp_name, int * send_comp_name_len,
				   char * send_grid_name, int * send_grid_name_len,
				   char * recv_comp_name, int * recv_comp_name_len,
				   char * recv_grid_name, int * recv_grid_name_len,
				   int * map_tag,
				   int * send_index, int * recv_index, double * coef, int * ngrid) ; 

void h3oup_set_interpolation_table_no_index(char * my_name, int * my_name_len,
					    char * send_comp_name, int * send_comp_name_len,
					    char * send_grid_name, int * send_grid_name_len,
					    char * recv_comp_name, int * recv_comp_name_len,
					    char * recv_grid_name, int * recv_grid_name_len,
					    int * map_tag) ; 
void h3oup_init_time(int * time_array, int * array_size) ;

void h3oup_set_time(char * my_name, int * my_name_len, int * time_array, int * array_size, int * delta_t) ;

void h3oup_get_vlayer(char * data_name, int * name_len, int * vlayer) ;

void h3oup_put_data_1d(char * data_name, int * name_len, double * data, int * data_len) ; 

void h3oup_put_data_25d(char * data_name, int * name_len, double * data, int * data_len1, int * data_len2) ;

void h3oup_get_data_1d(char * data_name, int * name_len, double * data, int * data_len, int * is_get_ok) ;

void h3oup_get_data_25d(char * data_name, int * name_len, double * data, int * data_len1, int * data_len2, int * is_get_ok) ;

void h3oup_get_num_of_put_data(int * num_of_data) ;

void h3oup_get_put_data_name(int * num_of_data, char * name, int * name_len) ;

void h3oup_get_num_of_get_data(int * num_of_data) ;

void h3oup_get_get_data_name(int * num_of_data, char * name, int * name_len) ;

void h3oup_send_array_int(char * my_comp_name, int * my_name_len, char * recv_comp_name, int * recv_name_len, int * array, int * array_size) ; 

void h3oup_recv_array_int(char * my_comp_name, int * my_name_len, char * send_comp_name, int * send_name_len, int * array, int * array_size) ; 

void h3oup_send_array_double(char * my_comp_name, int * my_name_len, char * recv_comp_name, int * recv_name_len, double * array, int * array_size) ; 

void h3oup_recv_array_double(char * my_comp_name, int * my_name_len, char * send_comp_name, int * send_name_len, double * array, int * array_size) ; 

void h3oup_coupling_end(int * time_array, int * nsize) ;

void h3oup_inc_calendar(int * time_array, int * nsize, int * delta_t)  ; 

void h3oup_init_simple(char * comp_name, int * name_len, char * log_level, int * log_level_len, int * stop_step, int * debug_mode) ;

void h3oup_end(void) ;

void h3oup_send_int_scalar(char * target_name, int * name_len, int * val) ;

void h3oup_send_real_scalar(char * target_name, int * name_len, float * val) ;

void h3oup_send_double_scalar(char * target_name, int * name_len, double * val) ;

void h3oup_send_int_array(char * target_name, int * name_len, int * val, int *array_size) ;

void h3oup_send_real_array(char * target_name, int * name_len, float * val, int * array_size) ;

void h3oup_send_double_array(char * target_name, int * name_len, double * val, int * array_size) ;

void h3oup_recv_int_scalar(char * source_name, int * name_len, int * val) ; 

void h3oup_recv_real_scalar(char * source_name, int * name_len, float * val) ; 

void h3oup_recv_double_scalar(char * source_name, int * name_len, double * val) ; 

void h3oup_recv_int_array(char * source_name, int * name_len, int * val, int * array_size) ; 

void h3oup_recv_real_array(char * source_name, int * name_len, float * val, int * array_size) ; 

void h3oup_recv_double_array(char * source_name, int * name_len, double * val, int * array_size) ; 
