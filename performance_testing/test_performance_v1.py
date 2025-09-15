import time
from v1_change_of_variable import festim_sim_v1_COV
import numpy as np

times = []

for i in range(10):
    start_time = time.perf_counter()
    festim_sim_v1_COV(n=50)
    end_time = time.perf_counter()
    times.append(end_time - start_time)

times = times[1:]
avg_time = np.mean(times)
print(f"Average time Change of Variable: {avg_time:.3f} s")

# save to a csv file
data = np.array(
    [
        ["Method", "Average Time (s)"],
        ["Change of Variable", f"{avg_time:.2f}"],
    ]
)
fname = "results/performance_comparison_v1.csv"
np.savetxt(fname, data, delimiter=",", fmt="%s")
print(f"Results saved to {fname}")
