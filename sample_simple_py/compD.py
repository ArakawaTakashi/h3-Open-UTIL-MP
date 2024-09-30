import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), "/work/jh210022a/share/h3-Open-UTIL-MP/h3ou/include/aquarius"))
import h3ou as h3opp
import numpy as np


my_name = "compD"

# initialize
h3opp.h3ou_init(my_name, "coupling.conf")

# get mpi parameter

my_size = h3opp.h3ou_get_my_size()
my_rank = h3opp.h3ou_get_my_rank()

print("compD my size = ", my_size, ", my_rank = ", my_rank)

# set my array
array_size = 10

my_array_int    = np.arange(1, array_size + 1, dtype = np.int32)
my_array_real   = np.arange(1, array_size + 1, dtype = np.float32)
my_array_double = np.arange(1, array_size + 1, dtype = np.float64)
target_array_int    = np.zeros(array_size, dtype = np.int32)
target_array_real   = np.zeros(array_size, dtype = np.float32)
target_array_double = np.zeros(array_size, dtype = np.float64)

# bcast global

target_array = h3opp.h3ou_bcast_global_int("compA", target_array_int)

print("bcast from compA = ", target_array)

target_array = h3opp.h3ou_bcast_global_real("compB", target_array_real)

print("bcast from compB = ", target_array)

target_array = h3opp.h3ou_bcast_global_double("compD", my_array_double)

print("bcast from compD = ", target_array)

# bcast model

target_array = h3opp.h3ou_bcast_model_real("compB", "compD", target_array_real)

print("compD bcast model from compC = ", target_array)

target_array = h3opp.h3ou_bcast_model_double("compD", "compB", my_array_double)

# send/recvt model

if (my_rank == 0):
    target_array = h3opp.h3ou_recv_model_real("compB", 0, target_array_real)
    print("compD send/recv model, ", target_array)
    h3opp.h3ou_send_model_double("compB", 0, my_array_double)

# bcast local

target_array = h3opp.h3ou_bcast_local_double(0, my_array_double)

print("compD bcast local, ", target_array)

# send/recv local

if (my_size >= 2):
    if (my_rank == 0):
        h3opp.h3ou_send_local_double(1, my_array_double)
    if (my_rank == 1):
        target_array = h3opp.h3ou_recv_local_double(0, target_array_double)
        print("compD send/recv local, ", target_array)
        
h3opp.h3ou_end()


