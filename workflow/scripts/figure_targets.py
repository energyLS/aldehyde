# Create graph of Moroccos CO2 emissions and climate targets

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import logging

logger = logging.getLogger(__name__)


def plot_figure_targets():
    """Plot graph showing historic emission and climate targets
    """
    url = 'https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv'
    emissions = pd.read_csv(url, index_col=0)

    cols = ["year","co2", "cement_co2", "total_ghg"]
    morocco_em = emissions[emissions.iso_code=='MAR'][cols].set_index("year")

    # Set up figure
    fig = plt.figure(facecolor='white')

    # Create Subplot TODO read emissions targets from config file
    ax = fig.add_subplot(1,1,1)
    ax.plot("co2", data=morocco_em, label='$\mathregular{CO_2}$ total', color="black")
    ax.plot("total_ghg", data=morocco_em, label='GHG total', color="grey")
    ax.scatter(2030, 115, label="2030 unconditional", marker = 'o', color = "green")
    ax.scatter(2030, 75, label="2030 conditional", marker = 'o', color = "green", alpha=0.7)

    ax.set_xlim(1940, 2060)
    ax.set_ylim(0, 120)

    # Further elements
    ax.set_ylabel("Emissions in Mt/year")
    ax.set_xlabel("Year")
    plt.legend()
    plt.grid(alpha=0.4)
    plt.tight_layout

    # Create graph and save it
    logger.info(f"save graph in {snakemake.output.figure}")
    plt.savefig(snakemake.output.figure)
    #plt.show()

    return

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake('figure_targets')

        sets_path_to_root('aldehyde')

    plot_figure_targets()