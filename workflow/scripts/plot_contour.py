import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging
import os
from matplotlib.colors import LinearSegmentedColormap

logger = logging.getLogger(__name__)


def prepare_data(data, zerofilter=False):
    # Prepare the data
    # Extract the Co2L which is in position 5-8
    data["opts"] = data["opts"].str[4:]

    # Rescale the cost from euro to B€
    data["cost"] = data["cost"] / 1e9

    # Round the data
    to_round = plottype
    data[to_round] = data[to_round].round(2)

    # Filter the data to remove 0 export and 0 co2 reduction
    print(f"zerofiler is set to {zerofilter} or in boolean {bool(zerofilter)}")
    if zerofilter:
        print("Filtering data")
        # data = data[(data["h2export"] != 0) & (data["opts"] != "2.0")]
        data = data[(data["h2export"] != 0)]

    else:
        pass

    return data


def reshape_data(data, opts, h2export, plottype):
    # Reshape the data for all columns in to_round and save it in data dictionary
    data_reshaped = {}
    # TODO loop is not necessary since I only use one plottype
    for i in [plottype]:
        data_reshaped[i] = data[i].values.reshape(len(opts), len(h2export)).T

    return data_reshaped


def plot_data(data_reshaped, plottype, levels, show_minimums, el_base_demand):
    # Turn "limit" to "reduction" (e.g. Co2L0.90 means 10% reduction)
    opts_reverse = 1 - opts
    opts_reverse[opts_reverse < 0] = 0

    # Plot a contour plot of the data having the y-axis the column "h2export", x-axis the column "sopts", and the z-axis the column "cost"
    fig = plt.figure(figsize=(9, 6))

    try:
        vmin = (snakemake.config["plot"]["contour_plot"]["vcontrol"][plottype][0],)
        vmax = (snakemake.config["plot"]["contour_plot"]["vcontrol"][plottype][1],)
    except:
        vmin = None
        vmax = None

    if snakemake.config["plot"]["contour_plot"]["normalize"] == False:
        if vmax is not None:
            contour = plt.contourf(
                opts_reverse * 100,
                h2export,
                np.flip(data_reshaped[plottype], axis=1),
                levels=levels,
                vmax=vmax[0],
                vmin=vmin[0],
            )
        else:
            contour = plt.contourf(
                opts_reverse * 100,
                h2export,
                np.flip(data_reshaped[plottype], axis=1),
                levels=levels,
            )
    elif snakemake.config["plot"]["contour_plot"]["normalize"] == True:
        # Option 1 : One end max and white at 0
        max_abs_value = np.max(np.abs(data_reshaped[plottype]))

        # Option 2: Both ends predefined and white at 0
        # max_abs_value = 60 # could also be obtained from config

        # (complement option 1 or 2)
        vmin = -max_abs_value
        vmax = max_abs_value
        cmap = plt.cm.RdYlGn_r

        # Option 3: Full range and white at 0
        # max_value = np.max(data_reshaped[plottype])
        # min_value = np.min(data_reshaped[plottype])
        # # norm = plt.Normalize(vmin=min_value, vmax=max_value)
        # zeropos = 1 - abs(max_value) / (abs(min_value) + abs(max_value))
        # cmap = LinearSegmentedColormap.from_list(
        #     "custom_cmap", [(0, "green"), (zeropos, "white"), (1, "red")]
        # )

        contour = plt.contourf(
            opts_reverse * 100,
            h2export,
            np.flip(data_reshaped[plottype], axis=1),
            levels=levels,
            vmax=vmax,
            vmin=vmin,
            cmap=cmap,
        )

    if show_minimums == True:
        # Return the position where the minimum of data_reshaped[plottype] is
        minpos = np.argmin(data_reshaped[plottype], axis=0)

        # Plot the minimum value as a black dot
        plt.plot(
            opts_reverse[::-1] * 100,
            h2export[minpos],
            "ko",
            markersize=4,
            label="min.",
            alpha=0.5,
        )

        # Plot approximation/regression with np.polyfit of the minimum value as line
        polydegree = 4
        line = np.poly1d(np.polyfit(opts_reverse * 100, h2export[minpos], polydegree))(
            np.linspace(min(opts_reverse), max(opts_reverse), 100) * 100
        )
        # Limit line to h2 export boundaries
        line[line < min(h2export)] = min(h2export)
        line[line > max(h2export)] = max(h2export)

        plt.plot(
            # opts_reverse[::-1] * 100,
            np.linspace(min(opts_reverse), max(opts_reverse), 100)[::-1] * 100,
            line,
            # "k--",
            color="black",
            label="min. approx.",
            alpha=0.5,
        )
        plt.legend()

    # plt.xlabel("CO$_2$ Reduction in % of base levels")
    plt.xlabel("Domestic mitigation in %")
    plt.ylabel("Hydrogen Export Volume in TWh")
    if snakemake.config["plot"]["contour_plot"]["normalize"] == False:
        plt.colorbar(pad=0.14).set_label(
            snakemake.config["plot"]["contour_plot"]["label"][plottype]
        )
    elif snakemake.config["plot"]["contour_plot"]["normalize"] == True:
        plt.colorbar(pad=0.14).set_label(
            snakemake.config["plot"]["contour_plot"]["norm_specs"][plottype]["label"]
        )
    else:
        pass

    ax2 = plt.gca().twinx()
    h2export_secondary = h2export / el_base_demand
    ax2.set_ylabel("H2-Exp/El-base-demand")
    ax2.set_ylim(0, max(h2export_secondary))  # Adjust the limit based on your data

    # Save the plot
    plt.savefig(snakemake.output.contour_plot, bbox_inches="tight")

    plt.show()

    # Set negative values in array df to 0
    # df[df < 0] = 0

    return


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "plot_contour",
            plottype="exp_AC_exclu_H2 El_all",
            levels=20,
            zerofilter="False",
        )

        sets_path_to_root("aldehyde")

    # Get the data
    data = pd.read_csv(snakemake.input.stats)

    # Get the plottype and levels (relevant for the contour plot)
    plottype = snakemake.wildcards.plottype
    zerofilter = snakemake.config["plot"]["contour_plot"]["zerofilter"]
    levels = int(snakemake.wildcards.levels)

    # Prepare the data
    data = prepare_data(data, zerofilter=zerofilter)

    # Get the unique values for opts and h2export required for the reshaping
    h2export = np.unique(data["h2export"])
    opts = np.unique(
        data["opts"].fillna(100).astype(float)
    )  # TODO improve the fillna value

    # Reshape the data
    data_reshaped = reshape_data(data, opts, h2export, plottype=plottype)

    # Plot the data
    show_minimums = snakemake.config["plot"]["contour_plot"]["show_minimums"]
    el_base_demand = min(
        reshape_data(data, opts, h2export, plottype="el_base_demand")["el_base_demand"][
            0
        ]
    )
    logger.info(f"el_base_demand is {el_base_demand} for all runs")

    if snakemake.config["plot"]["contour_plot"]["normalize"] == True:
        if (
            snakemake.config["plot"]["contour_plot"]["norm_specs"][plottype][
                "normalize_by"
            ]
            == "export"
        ):
            data_reshaped[plottype] = (
                (data_reshaped[plottype] / data_reshaped[plottype][0, :][np.newaxis, :])
                - 1
            ) * 100
        elif (
            snakemake.config["plot"]["contour_plot"]["norm_specs"][plottype][
                "normalize_by"
            ]
            == "decarb"
        ):
            data_reshaped[plottype] = (
                (data_reshaped[plottype] / data_reshaped[plottype][:, 0][:, np.newaxis])
                - 1
            ) * 100
        elif (
            snakemake.config["plot"]["contour_plot"]["norm_specs"][plottype][
                "normalize_by"
            ]
            == "exdecarb"
        ):
            data_reshaped[plottype] = (
                (data_reshaped[plottype] / data_reshaped[plottype][0, 0]) - 1
            ) * 100
        else:
            pass

    else:
        pass

    plot_data(
        data_reshaped,
        plottype,
        levels,
        show_minimums=show_minimums,
        el_base_demand=el_base_demand,
    )
