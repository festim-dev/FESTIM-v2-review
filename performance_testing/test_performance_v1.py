import time
from v1_change_of_variable import festim_sim_v1_COV
import numpy as np

mesh_sizes = np.geomspace(10, 1000, num=4, dtype=int)
print(mesh_sizes)
exit()
times = []

for i in range(10):
    start_time = time.time()
    festim_sim_v1_COV(n=100)
    end_time = time.time()
    times.append(end_time - start_time)

times = times[1:]
avg_time = np.mean(times)
print(f"Average simulation time: {avg_time:.3f} s")
