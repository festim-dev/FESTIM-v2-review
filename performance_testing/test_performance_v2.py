import time
from v2_change_of_variable import festim_sim_v2_COV
from v2_discontinuous_penalty import festim_sim_v2_disc_pen
from v2_discontinuous_nitsche import festim_sim_v2_disc_nitsche
import numpy as np

from dolfinx.log import set_log_level, LogLevel

# set_log_level(LogLevel.INFO)
mesh_size = 50
N = 10
print("running penalty cases...")
times = []
for i in range(N):
    start_time = time.perf_counter()
    festim_sim_v2_disc_pen(n=mesh_size)
    end_time = time.perf_counter()
    times.append(end_time - start_time)

avg_time_2 = np.mean(times[1:])

print("running nitsche cases...")
times = []
for i in range(N):
    start_time = time.perf_counter()
    festim_sim_v2_disc_nitsche(n=mesh_size)
    end_time = time.perf_counter()
    times.append(end_time - start_time)

avg_time_3 = np.mean(times[1:])

print("running change of variable cases...")
times = []
for i in range(N):
    start_time = time.perf_counter()
    festim_sim_v2_COV(n=mesh_size)
    end_time = time.perf_counter()
    times.append(end_time - start_time)

avg_time_1 = np.mean(times[1:])

print(f"Average time Discontinuous Penalty: {avg_time_2:.2f} seconds")
print(f"Average time Discontinuous Nitsche: {avg_time_3:.2f} seconds")
print(f"Average time Change of Variable: {avg_time_1:.2f} seconds")

# save to a csv file
data = np.array(
    [
        ["Method", "Average Time (s)"],
        ["Change of Variable", f"{avg_time_1:.2f}"],
        ["Discontinuous Penalty", f"{avg_time_2:.2f}"],
        ["Discontinuous Nitsche", f"{avg_time_3:.2f}"],
    ]
)
fname = "results/performance_comparison_v2.csv"
np.savetxt(fname, data, delimiter=",", fmt="%s")
print(f"Results saved to {fname}")
