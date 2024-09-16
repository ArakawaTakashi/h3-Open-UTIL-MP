import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), "/home/arakawa/work/2024H3OPEN/h3ou2024/h3-Open-UTIL-MP/h3ou/include"))
import h3ou as h3opp
import numpy as np

h3opp.h3ou_init("compD", "coupling.conf")

my_rank = h3opp.h3ou_get_my_rank()
print("my_rank = ", my_rank)

my_size = h3opp.h3ou_get_my_size()
print("my_size = ", my_size)


target_data = 0
my_data = 100000
recv_data = h3opp.h3ou_recv_scalar("compA", target_data)
h3opp.h3ou_send_scalar("compA", my_data)

print("compD target_data = ", recv_data)

send_array = [1, 2, 3, 4, 5]
recv_array = [0, 0, 0, 0, 0]

h3opp.h3ou_send_array("compC", send_array)
recv_array = h3opp.h3ou_recv_array("compC", recv_array)


print("compD recv_array = ", recv_array)

recv_array = [0, 0, 0, 0, 0]
 
irecv_array = np.zeros(10, dtype = np.int32)

source_pe = 0

if (my_rank ==1):
    recv_list = h3opp.h3ou_recv_model_int("compA", source_pe, irecv_array)
    irecv_array[:] = recv_list[:]
    print("irecv call OK ", recv_list)

recv_list = h3opp.h3ou_bcast_local_int(1, irecv_array)
print("bcast call OK ", recv_list)


if (my_rank == 0):
    h3opp.h3ou_send_local_int(1, send_array)
elif (my_rank == 1):
    recv_array = h3opp.h3ou_recv_local_int(0, recv_array)
    print("h3ou_recv_local, recv_array = ", recv_array)
    
h3opp.h3ou_end()


