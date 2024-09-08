import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), "/home/arakawa/work/2024H3OPEN/h3ou2024/h3ou/include"))
import h3ou as h3opp
import numpy as np

h3opp.h3ou_init_simple("compD", "SILENT", 0, 0)

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


h3opp.h3ou_end()


