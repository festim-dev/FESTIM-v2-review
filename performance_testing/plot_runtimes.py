import numpy as np
import matplotlib
from pathlib import Path
matplotlib.use("Agg")  # Use non-interactive backend
import matplotlib.pyplot as plt

folder = Path("results")

v1_data = np.genfromtxt(
    folder/"performance_comparison_v1.csv", delimiter=",", names=True
)
v2_data = np.genfromtxt(
    folder/"performance_comparison_v2.csv", delimiter=",", names=True
)

v1_cov_data = float(v1_data["Average_Time_s"])

# Access v2 data by index since we know the order from the CSV
v2_cov_data = float(v2_data["Average_Time_s"][0])  # Change of Variable
v2_pen_data = float(v2_data["Average_Time_s"][1])  # Discontinuous Penalty
v2_niet_data = float(v2_data["Average_Time_s"][2])  # Discontinuous Nitsches

# Create bar chart with two groups: v1 (single bar) and v2 (three bars)
fig, ax = plt.subplots(figsize=(10, 6))

# Define positions for the bars
v1_pos = [0]  # Single position for v1
v2_pos = [2, 2.8, 3.6]  # Three positions for v2 methods

# V1 data (single bar)
v1_bar = ax.bar(
    v1_pos, [v1_cov_data], width=0.6, label="FESTIM v1", color="#D8D8D2", alpha=0.8
)

# V2 data (three bars)
v2_bars = ax.bar(
    v2_pos,
    [v2_cov_data, v2_pen_data, v2_niet_data],
    width=0.6,
    color=["#3A3F3E", "#2A5350", "#F7B614"],
    alpha=0.8,
)

# Customize the plot
# Remove y-label and y-ticks for cleaner look
ax.set_ylabel("")
ax.set_yticks([])

# Set x-axis labels for each bar
all_positions = v1_pos + v2_pos
all_labels = ["Change of\nVariable"] + [
    "Change of\nVariable",
    "Discontinuous\nPenalty",
    "Discontinuous\nNitsches",

]
ax.set_xticks(all_positions)
ax.set_xticklabels(all_labels, fontsize=9)

# Add FESTIM v1 and v2 annotations
ax.annotate(
    "FESTIM v1",
    xy=(0, v1_cov_data / 2),
    xytext=(0, v1_cov_data / 2),
    ha="center",
    va="center",
    fontsize=12,
    fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
)

ax.annotate(
    "FESTIM v2",
    xy=(2.8, max(v2_cov_data, v2_pen_data, v2_niet_data) * 1.5),
    xytext=(2.8, max(v2_cov_data, v2_pen_data, v2_niet_data) * 1.5),
    ha="center",
    va="center",
    fontsize=12,
    fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
)

# Add value labels on top of bars
ax.text(
    v1_pos[0],
    v1_cov_data + max(v2_cov_data, v2_pen_data, v2_niet_data) * 0.02,
    f"{v1_cov_data:.1f} s",
    ha="center",
    va="bottom",
)

for pos, val in zip(v2_pos, [v2_cov_data, v2_pen_data, v2_niet_data]):
    ax.text(
        pos,
        val + max(v2_cov_data, v2_pen_data, v2_niet_data) * 0.02,
        f"{val:.1f} s",
        ha="center",
        va="bottom",
    )

# Remove top, right, and left spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

# Adjust layout and save
plt.tight_layout()
plt.savefig(
    "results/performance_comparison_bar_chart.pdf", dpi=300, bbox_inches="tight"
)
# print("Bar chart saved as 'results/performance_comparison_bar_chart.png'")
# plt.show()  # Commented out to avoid display issues
