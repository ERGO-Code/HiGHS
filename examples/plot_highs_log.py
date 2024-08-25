import re
import matplotlib.pyplot as plt
import numpy as np


def parse_highs_log(log_file_path):
    last_full_entry = []
    current_entry = []
    found_solution = False

    with open(log_file_path, "r") as f:
        for line in f:
            if "Running HiGHS" in line:
                if found_solution:
                    last_full_entry = current_entry
                current_entry = [line]
                found_solution = False
            else:
                current_entry.append(line)
                if "Writing the solution to" in line:
                    found_solution = True

    if not last_full_entry:
        last_full_entry = current_entry

    if not last_full_entry:
        return None, None, None, None, None, None

    time_values, best_bound_values, best_sol_values, in_queue_values, expl_values, gap_values = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for line in last_full_entry:
        match = re.search(r"\dk?\s+\d+\.\ds$", line)

        if not match:
            continue

        tokens = line.split()
        if len(tokens) == 13:
            tokens = tokens[1:]
        assert len(tokens) == 12, f"{line}"

        in_queue_values.append(float(tokens[1]))  # InQueue
        expl_values.append(float(tokens[3].replace("%", "")))  # Expl.%
        best_bound_values.append(float(tokens[4].replace("inf", "nan")))  # Best Bound
        best_sol_values.append(float(tokens[5].replace("inf", "nan")))  # Best Sol
        gap_values.append(
            float(tokens[6].replace("%", "").replace("inf", "nan").replace("Large", "nan"))
        )  # Gap%
        time_values.append(float(tokens[11].replace("s", "")))  # Time

    return time_values, best_bound_values, best_sol_values, in_queue_values, expl_values, gap_values


def plot_highs_log(
    time_values, best_bound_values, best_sol_values, in_queue_values, expl_values, gap_values
):
    fig, ax1 = plt.subplots(figsize=(10, 6))

    best_bound_colour = "blue"
    best_solution_colour = "green"
    in_queue_colour = "red"
    explored_colour = "purple"
    gap_colour = "orange"

    # Plot Objective Bounds
    ax1.plot(time_values, best_bound_values, label="Best Bound", color=best_bound_colour)
    ax1.plot(time_values, best_sol_values, label="Best Solution", color=best_solution_colour)
    ax1.set_xlabel("Time (seconds)")
    ax1.set_ylabel("Objective Bounds", color=best_bound_colour, labelpad=15)
    ax1.tick_params(axis="y", labelcolor=best_bound_colour)

    # Limit y-axis to the range between min and max of the non-NaN values
    valid_gap_index = next(i for i, gap in enumerate(gap_values) if not np.isnan(gap))
    min_y = min(best_bound_values[valid_gap_index], best_sol_values[valid_gap_index])
    max_y = max(best_bound_values[valid_gap_index], best_sol_values[valid_gap_index])
    padding = (max_y - min_y) * 0.1
    ax1.set_ylim(min_y - padding, max_y + padding)

    # Add second y-axis for InQueue values
    ax2 = ax1.twinx()
    ax2.plot(time_values, in_queue_values, label="InQueue", color=in_queue_colour)
#    ax2.set_ylabel("InQueue", color=in_queue_colour, loc="top", labelpad=12)
    ax2.yaxis.label.set_rotation(0)
    ax2.tick_params(axis="y", labelcolor=in_queue_colour)

    # Add third y-axis for Explored % values (scaled)
    ax3 = ax1.twinx()
    ax3.spines["right"].set_position(("outward", 50))
    ax3.plot(time_values, expl_values, label="Explored %", color=explored_colour)
#    ax3.set_ylabel("Expl.%", color=explored_colour, loc="top", labelpad=10)
    ax3.yaxis.label.set_rotation(0)
    ax3.tick_params(axis="y", labelcolor=explored_colour)

    # Add fourth y-axis for Gap % values (scaled)
    ax4 = ax1.twinx()
    ax4.spines["right"].set_position(("outward", 90))
    ax4.plot(time_values, gap_values, label="Gap %", color=gap_colour, linestyle="--", linewidth=0.5)
#    ax4.set_ylabel("Gap.%", color=gap_colour, loc="top", labelpad=22)
    ax4.yaxis.label.set_rotation(0)
    ax4.tick_params(axis="y", labelcolor=gap_colour)

    # Plot vertical hash lines where Best Solution changes
    for i in range(1, len(best_sol_values)):
        if best_sol_values[i] != best_sol_values[i - 1]:  # Change detected
            ax1.axvline(x=time_values[i], color="grey", linestyle="--", linewidth=0.5)

    # Shift plot area left to make room on the right for the three y-axis labels.
    fig.subplots_adjust(left=0.08, right=0.85)

    # Set up legend
    fig.legend(loc="lower center", ncols=5)

    # Show plot
    plt.title("HiGHS MIP Log Analysis", pad=20)
    plt.show()


#log_file_path = "/path/to/your/logfile.log"
log_file_path = "HiGHS.log"
time_values, best_bound_values, best_sol_values, in_queue_values, expl_values, gap_values = (
    parse_highs_log(log_file_path)
)

plot_highs_log(
    time_values, best_bound_values, best_sol_values, in_queue_values, expl_values, gap_values
)
