# Create graph of Moroccos CO2 emissions and climate targets

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import os
import logging

logger = logging.getLogger(__name__)


def plot_figure_land_conflicts():

    morocco_gadm = gpd.read_file(snakemake.input.gadm_shapes)

    morocco_gadm["islanded_installation"] = (0.8 - morocco_gadm["pop"]*(1/max(morocco_gadm["pop"])) ) *100
    morocco_gadm["integrated_installation"] = (0.8 - morocco_gadm["gdp"]*(1/max(morocco_gadm["gdp"])) ) *100
    morocco_gadm["both_installation"] = (morocco_gadm["islanded_installation"] + morocco_gadm["integrated_installation"] )

    fig = plt.figure(facecolor='white', figsize=(10,7))
    fig.suptitle("Land conflicts")

    ax = fig.add_subplot(2,2,1, title="MAR ES only")
    ax2 = fig.add_subplot(2,2,3, title="Export only")
    ax3 = fig.add_subplot(2,2,2, title="MAR ES + Export")

    l1=morocco_gadm.plot(column='islanded_installation', ax=ax, cmap='OrRd', legend=True, vmax=200, legend_kwds={'label': "Land use in %"})  #, scheme='quantiles'
    l2=morocco_gadm.plot(column='integrated_installation', ax=ax2, cmap='OrRd', legend=True, vmax=200, legend_kwds={'label': "Land use in %"})  #, scheme='quantiles'
    l3=morocco_gadm.plot(column='both_installation', ax=ax3, legend=True, cmap='OrRd', vmin=0, vmax=200, legend_kwds={'label': "Land use in %"})  #, scheme='quantiles'

    ax.set_axis_off()
    ax2.set_axis_off()
    ax3.set_axis_off()

    # Create graph and save it
    logger.info(f"save graph in {snakemake.output.land_conflicts}")
    plt.savefig(snakemake.output.land_conflicts)

    return

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
    ax.scatter(2030, snakemake.config["climate_targets"]["2030_cond"], label="2030 unconditional", marker = 'o', color = "green")
    ax.scatter(2030, snakemake.config["climate_targets"]["2030_uncond"], label="2030 conditional", marker = 'o', color = "green", alpha=0.7)

    ax.set_xlim(1940, 2100)
    ax.set_ylim(0, 120)

    # Further elements
    ax.set_ylabel("Emissions in Mt/year")
    ax.set_xlabel("Year")
    plt.legend()
    plt.grid(alpha=0.4)
    plt.tight_layout

    # Create graph and save it
    logger.info(f"save graph in {snakemake.output.morocco_em}")
    plt.savefig(snakemake.output.morocco_em)

    return

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake('plot_figures')

        sets_path_to_root('aldehyde')

    plot_figure_land_conflicts()

    plot_figure_targets()