# Create graph of Moroccos CO2 emissions and climate targets

import pypsa
import yaml
import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import cartopy.crs as ccrs
import cartopy

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from matplotlib.patches import Circle, Patch
from matplotlib.legend_handler import HandlerPatch

from pypsa.descriptors import get_switchable_as_dense as as_dense

from shapely import wkt
import sys, os

import logging

logger = logging.getLogger(__name__)

def plot_spatial(
    type,
    df,
    geodf,
    carrier,
    cmap="Blues",
    vmax=100,
    vmin=0,
    #label="capacity factors [%]",
    fn=None,
):
    label_type = label[type]

    proj = ccrs.EqualEarth()
    geodf = geodf.to_crs(proj.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": proj})

    geodf.plot(
        ax=ax,
        column=df[carrier].reindex(geodf.index),
        cmap=cmap,
        linewidths=0,
        legend=True,
        vmax=vmax,
        vmin=vmin,
        legend_kwds={
            "label": label_type,
            "shrink": 0.7,
            # "extend": "max",
        },
    )

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.2, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.2, zorder=2)

    # Add title
    ax.set_title(f"{snakemake.wildcards.network} (Export {snakemake.wildcards.h2export} TWh)")

    #plt.gca().outline_patch.set_visible(False)
    ax.set_facecolor("white")

    print(snakemake.output.plot)
    plt.savefig(snakemake.output.plot, bbox_inches="tight")


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake('plot_spatial_figures', energy_source = "solar", plottype = "p-nom-max", h2export=20)

                    

        sets_path_to_root('aldehyde')

    label = {"cf": "Capacity factors in %",
        "p-nom-opt": "Optimised capacity in GW",
        "p-nom-max": "Potentials in GW",
        "land-use": "Land use of available potential in %",
        }

    nodes = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    # Import the network 
    #n = pypsa.Network(snakemake.input.network)
    n = pypsa.Network(snakemake.input["network_"+ snakemake.wildcards.network])

    # Insert hack to copy index of all rows of n.buses with carrier "AC" to the column "location" and avoid length mismatch. They were missing
    n.buses.loc[n.buses.carrier == "AC", "location"] = n.buses.index[n.buses.carrier == "AC"]

    techs = [
        "onwind",
        "solar",
    ]

    energy_source = snakemake.wildcards.energy_source
    plottype = snakemake.wildcards.plottype

    if plottype == "cf":
        # Create a dataframe with the mean capacity factors of the generators
        df = (
            n.generators_t.p_max_pu.mean()
            .groupby([n.generators.carrier, n.generators.bus.map(n.buses.location)])
            .first()
            .unstack(0)
            .mul(100)
            )
        df = df[techs]

        if energy_source == "onwind":
            plot_spatial(plottype, df, nodes, energy_source, cmap="Blues", vmax=55, vmin=0)

        elif energy_source == "solar":
            plot_spatial(plottype, df, nodes, energy_source, cmap="Oranges", vmax=max(df["solar"]), vmin=0)


    elif plottype == "p-nom-opt":

        df = (
            n.generators.p_nom_opt
            .groupby([n.generators.carrier, n.generators.bus.map(n.buses.location)])
            .first()
            .unstack(0)
            .div(1e3) # in GW
            )
        if energy_source == "onwind":
            plot_spatial(plottype, df, nodes, "onwind", vmax=max(df["onwind"]))
        elif energy_source == "solar":
            plot_spatial(plottype, df, nodes, "solar", cmap="Oranges", vmax=max(df["solar"]), vmin=0)

    elif plottype == "p-nom-max":
        df = (
            n.generators.p_nom_max
            .groupby([n.generators.carrier, n.generators.bus.map(n.buses.location)])
            .first()
            .unstack(0)
            .div(1e3) # in GW
            )
        if energy_source == "onwind":
            plot_spatial(plottype, df, nodes, energy_source, vmax=max(df[energy_source]))
        elif energy_source == "solar":
            plot_spatial(plottype, df, nodes, energy_source, cmap="Oranges", vmax=max(df[energy_source]), vmin=0)

    elif plottype == "land-use":
        df = (
            (n.generators.p_nom_opt/n.generators.p_nom_max)
            .groupby([n.generators.carrier, n.generators.bus.map(n.buses.location)])
            .first()
            .unstack(0)
            .mul(1e2) # in %
            )
        if energy_source == "onwind":
            plot_spatial(plottype, df, nodes, energy_source, vmax=100)
        elif energy_source == "solar":
            plot_spatial(plottype, df, nodes, energy_source, cmap="Oranges", vmax=100)
