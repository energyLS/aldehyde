import pypsa
import os
import pandas as pd
import geopandas as gpd
from _helpers import load_costs
import logging

logger = logging.getLogger(__name__)

def add_hydrogen(n, costs):
    "function to add hydrogen as an energy carrier with its conversion technologies from and to AC"

    n.add("Carrier", "H2")

    n.madd(
        "Bus",
        nodes + " H2",
        location=nodes,
        carrier="H2",
        x=n.buses.loc[list(nodes)].x.values,
        y=n.buses.loc[list(nodes)].y.values,
    )

    n.madd(
        "Link",
        nodes + " H2 Electrolysis",
        bus1=nodes + " H2",
        bus0=nodes,
        p_nom_extendable=True,
        carrier="H2 Electrolysis",
        efficiency=costs.at["electrolysis", "efficiency"],
        capital_cost=costs.at["electrolysis", "capital_cost"],
        lifetime=costs.at["electrolysis", "lifetime"],
    )

def add_export(n, h2export):
    "function to add export of hydrogen to the network"

    # Read country shape
    country_shape = gpd.read_file(snakemake.input["shapes_path"])
    # Find most northwestern point in country shape and get x and y coordinates
    country_shape = country_shape.to_crs("EPSG:4326")

    # Get coordinates of the most western and northern point of the country and add a buffer of 2 degrees (equiv. to approx 220 km)
    x_export = country_shape.geometry.centroid.x.min() -2
    y_export = country_shape.geometry.centroid.y.max() +2

    # add export bus
    n.add(
        "Bus",
        "H2 export bus",
        carrier="H2",
        x=x_export,
        y=y_export,
    )

    # add export links
    logger.info("Adding export links")
    n.madd(
        "Link",
        names=nodes + " H2 export",
        bus0=nodes + " H2",
        bus1="H2 export bus",
        p_nom_extendable=True,
    )

    export_links = n.links[n.links.index.str.contains("export")]
    logger.info(export_links)

    # add store
    n.add(
        "Store",
        "H2 export store",
        bus="H2 export bus",
        e_nom_extendable=True,
        carrier="H2",
        e_initial=0,
        marginal_cost=0,
        capital_cost=0,
    )

    # add load
    n.add(
        "Load",
        "H2 export load",
        bus="H2 export bus",
        carrier="H2",
        p_set=h2export*1e6 / 8760, # convert TWh to MWh and divide by 8760 to get hourly average
    )

    return

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "add_esc",
            simpl="",
            clusters="4",
            ll="c1.0",
            opts="Co2L0.10",
            planning_horizons="2030",
            sopts="144H",
            discountrate=0.071,
            demand="DF",
            h2export=20,
        )

        sets_path_to_root('aldehyde')


    # Read Morocco network
    # https://pypsa.readthedocs.io/en/latest/components.html?highlight=override_component_attrs#custom-components
    n = pypsa.Network(snakemake.input.network) #, override_component_attrs=overrides)

    Nyears = n.snapshot_weightings.generators.sum() / 8760

    costs = load_costs(
        snakemake.input.tech_costs,
        snakemake.config["costs"],
        snakemake.config["electricity"],
        Nyears,
    )

    #nodes = n.buses.index

    #nodes = n.buses[n.buses.carrier == "AC"].index  # TODO if you take nodes from the index of buses of n it's more than pop_layout
    # clustering of regions must be double checked.. refer to regions onshore


    # add hydrogen buses
    #add_hydrogen(n, costs)

    # get export demand from wildcard
    #h2export = eval(snakemake.wildcards["h2export"]) # in TWh


    # logger.info(
    #     f"The yearly export demand is {h2export} TWh resulting in an hourly average of {h2export*1e6/8760:.2f} MWh"
    # )

    # add export value and components to network
    #add_export(n, h2export)
    
    n.export_to_netcdf(snakemake.output[0])