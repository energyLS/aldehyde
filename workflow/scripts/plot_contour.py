import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging
import os

logger = logging.getLogger(__name__)


def prepare_data(data, zerofilter=False):

    # Prepare the data
    # Extract the Co2L which is in position 5-8
    data["opts"] = data["opts"].str[4:]

    #Rescale the cost from euro to B€
    data["cost"] = (data["cost"]/1e9)

    # Round the data
    to_round = plottype
    data[to_round] = data[to_round].round(2)

    # Filter the data to remove 0 export and 0 co2 reduction
    print(f"zerofiler is set to {zerofilter} or in boolean {bool(zerofilter)}")
    if zerofilter:
       print("Filtering data")
       data = data[(data["h2export"] != 0) & (data["opts"] != "2.0")]
    else:
        pass

    return data

def reshape_data(data, opts, h2export):
    # Reshape the data for all columns in to_round and save it in data dictionary
    data_reshaped = {}
    # TODO loop is not necessary since I only use one plottype
    for i in [plottype]:
        data_reshaped[i] = data[i].values.reshape(len(opts), len(h2export)).T

    return data_reshaped

def plot_data(data_reshaped, plottype,levels):

    # Turn "limit" to "reduction" (e.g. Co2L0.90 means 10% reduction)
    opts_reverse = 1-opts
    opts_reverse[opts_reverse < 0] = 0

    # Plot a contour plot of the data having the y-axis the column "h2export", x-axis the column "sopts", and the z-axis the column "cost"
    plt.contourf(opts_reverse*100,h2export,np.flip(data_reshaped[plottype], axis=1), levels=levels)
    plt.xlabel('CO$_2$ Reduction in % of base levels')
    plt.ylabel('Hydrogen Export Volume in TWh')
    plt.colorbar().set_label(snakemake.config["plot"]["contour_plot"]["label"][plottype])

    # Save the plot
    plt.savefig(snakemake.output.contour_plot, bbox_inches='tight')

    plt.show()

    # Set negative values in array df to 0
    # df[df < 0] = 0


    return

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake('plot_contour', plottype = "lcoh_w_mixed", levels=10, zerofilter=False)


        sets_path_to_root('aldehyde')

    # Get the data
    data = pd.read_csv(snakemake.input.stats)

    # Get the plottype and levels (relevant for the contour plot)
    plottype = snakemake.wildcards.plottype
    zerofilter = snakemake.config["plot"]["contour_plot"]["zerofilter"]
    levels = int(snakemake.wildcards.levels)

    # Prepare the data
    data = prepare_data(data, zerofilter=zerofilter)

    # Get the unique values for opts and h2export required for the reshaping
    h2export = np.unique(data['h2export'])
    opts = np.unique(data['opts'].fillna(100).astype(float)) #TODO improve the fillna value


    # Reshape the data
    data_reshaped = reshape_data(data, opts, h2export)

    # Plot the data
    plot_data(data_reshaped, plottype, levels)