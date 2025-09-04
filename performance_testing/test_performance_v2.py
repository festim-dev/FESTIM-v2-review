import time
from v2_change_of_variable import festim_sim_v2_COV
from v2_discontinuous_penalty import festim_sim_v2_disc_pen
from v2_discontinuous_nietsche import festim_sim_v2_disc_nietsche
import numpy as np

mesh_size = 500

print("running change of variable cases...")
times = []
for i in range(10):
    start_time = time.time()
    festim_sim_v2_COV(n=mesh_size)
    end_time = time.time()
    times.append(end_time - start_time)

avg_time_1 = np.mean(times[1:])

print("running penalty cases...")
times = []
for i in range(10):
    start_time = time.time()
    festim_sim_v2_disc_pen(n=mesh_size)
    end_time = time.time()
    times.append(end_time - start_time)

avg_time_2 = np.mean(times[1:])

print("running nietsche cases...")
times = []
for i in range(10):
    start_time = time.time()
    festim_sim_v2_disc_nietsche(n=mesh_size)
    end_time = time.time()
    times.append(end_time - start_time)

avg_time_3 = np.mean(times[1:])

print(f"Average time Change of Variable: {avg_time_1:.2f} seconds")
print(f"Average time Discontinuous Penalty: {avg_time_2:.2f} seconds")
print(f"Average time Discontinuous Nietsches: {avg_time_3:.2f} seconds")
