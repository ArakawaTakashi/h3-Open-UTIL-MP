import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), "/home/arakawa/work/2024H3OPEN/h3ou2024/h3-Open-UTIL-MP/h3ou/include"))
import h3ou as h3opp
import numpy as np
from datetime import datetime, timedelta


my_name = "compD"
global_array_size = 20
init_time = datetime(2000, 1, 1, 0, 0, 0)
delta_t = 3600

# initialize
h3opp.h3ou_init(my_name, "coupling.conf")

# get mpi parameter

my_size = h3opp.h3ou_get_my_size()
my_rank = h3opp.h3ou_get_my_rank()

print("compD my size = ", my_size, ", my_rank = ", my_rank)

if ((global_array_size % my_size) != 0):
    print("Error, mod(global_array_size, my_size) must be 0\n") ;
    sys.exit()

array_size = int(global_array_size/my_size)
print("compD array size = ", array_size)

my_grid = np.zeros((array_size), dtype=np.int32)

for i in range(array_size):
    my_grid[i] = my_size * my_rank + i + 1

h3opp.h3ou_def_grid(my_grid, my_name, "compD_grid1", 20)
h3opp.h3ou_end_grid_def()

h3opp.h3ou_set_interpolation_table_no_index(my_name, "compA", "compA_grid1", "compD", "compD_grid1", 1)
h3opp.h3ou_set_interpolation_table_no_index(my_name, "compD", "compD_grid1", "compA", "compA_grid1", 1)

time_array = np.array([init_time.year, init_time.month, init_time.day, init_time.hour, init_time.minute, init_time.second], dtype = np.int32)

h3opp.h3ou_init_time(time_array)


my_array_int    = np.arange(1, array_size + 1, dtype = np.int32)
my_array_real   = np.arange(1, array_size + 1, dtype = np.float32)
my_array_double = np.arange(1, array_size + 1, dtype = np.float64)
target_array_int    = np.zeros(array_size, dtype = np.int32)
target_array_real   = np.zeros(array_size, dtype = np.float32)
target_array_double = np.zeros(array_size, dtype = np.float64)


h3opp.h3ou_put_data_1d("compD_var1", my_array_double)

now_time = init_time

for t in range(24):
    time_array = np.array([now_time.year, now_time.month, now_time.day, now_time.hour, now_time.minute, now_time.second], dtype = np.int32)
    h3opp.h3ou_set_time(my_name, time_array, delta_t)

    is_ok = h3opp.h3ou_get_data_1d("compA_var3", target_array_double)

    if (is_ok):
        print("get compA_var3 ", is_ok, target_array_double)

    target_array_double[:] = 0.0
    
    h3opp.h3ou_put_data_1d("compD_var1", my_array_double)
    

    now_tome = now_time + timedelta(seconds = delta_t)

time_array = np.array([now_time.year, now_time.month, now_time.day, now_time.hour, now_time.minute, now_time.second], dtype=np.int32)
h3opp.h3ou_coupling_end(time_array)


