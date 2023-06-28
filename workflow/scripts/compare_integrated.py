"""Script to compare networks.
"""

import pypsa
import pandas as pd
from matplotlib.pyplot import legend
import matplotlib.pyplot as plt
import numpy as np
import os
import logging


def get_bus_demand(n, busname):
    """Get the demand at a certain bus (includes stores/StorageUnits) based on the energy_balance() function

    Parameters
    ----------
    busname : _type_
        _description_
    n : _type_
        _description_
    """
    energy_balance = n.statistics.energy_balance(
        aggregate_bus=False, aggregate_time=False
    )
    # Get the energy balance of the bus specified in busname
    energy_balance_bus = energy_balance.loc[:, :, :, busname]

    # Filter for negative values and sum them up (note: may include stores/StorageUnits)
    demand = energy_balance_bus[energy_balance_bus < 0]
    demand = demand.groupby("bus").sum()

    return demand


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("compare_integrated")

        sets_path_to_root("aldehyde")

    # Load networks from snakefile input
    # Create loop over snakemake inputs to load all networks

    # limit the number of networks
    # limit = 5

    # Loop over the networks and save the objective of the networks in array according to the h2export and opts wildcards
    cost_df = pd.DataFrame(columns=["h2export", "opts", "cost", "lcoh"])
    cost0_df = pd.DataFrame(columns=["h2export", "opts", "cost"])

    for i in snakemake.input.networks:  # [0:limit]:
        n = pypsa.Network(i)
        cost = n.objective

        # Get the h2export (before "export") and opts values from the network name
        h2export = i.split("_")[-1][:-9]
        opts = i.split("_")[-6]

        # Read and save the objective of the network if h2export is 0. Save the costs in a new array
        if h2export == "0":
            n0 = pypsa.Network(i)
            cost0 = n0.objective
            # Save the cost0 in an array
            cost0_df = pd.concat(
                [
                    cost0_df,
                    pd.DataFrame(
                        {"h2export": [h2export], "opts": [opts], "cost": [cost0]}
                    ),
                ],
                ignore_index=True,
            )

        # Calculate the cost difference between the networks n and n0
        cost_difference = cost - cost0_df[cost0_df["opts"] == opts]["cost"]  # in €

        # Calculate the cost of the H2 export by taking the difference in system costs and dividing by the export demand
        lcoh = cost_difference / (
            n.loads_t.p.loc[:, "H2 export load"].sum()
            * n.snapshot_weightings.generators[0]
        )  # Cost difference between systems divided by the export demand

        # Calculate the marginal price of the buses: Warning: not weighted price
        lcoh_marginal = (
            n.buses_t.marginal_price.mean().groupby(n.buses.carrier).mean().loc["H2"]
        )  # in €/MWh (correct calc, checked with pypsa-earth-sec results)
        lcoh_marginal_export = n.buses_t.marginal_price.mean().loc[
            "H2 export bus"
        ]  # in €/MWh
        lcoe_marginal = (
            n.buses_t.marginal_price.mean().groupby(n.buses.carrier).mean().loc["AC"]
        )  # in €/MWh (correct calc, checked with pypsa-earth-sec results)

        # Calculate the weighted marginal price of the buses
        # TODO doublecheck, why the weighted marginal price reuqires the links only and no loads
        buses = n.buses.index[
            (n.buses.index.str[-2:] == "H2") | (n.buses.index == "H2 export bus")
        ]
        load = pd.DataFrame(index=n.snapshots, columns=buses, data=0.0)

        for tech in ["Sabatier", "H2 Fuel Cell"]:
            names = n.links.index[n.links.index.to_series().str[-len(tech) :] == tech]

            load += (
                n.links_t.p0[names].groupby(n.links.loc[names, "bus0"], axis=1).sum()
            )

        lcoh_marginal_weighted = (
            load * n.buses_t.marginal_price[buses]
        ).sum().sum() / load.sum().sum()

        ############################################################
        # New approach to calculate the weighted marginal price of the hydrogen buses: sum(Marginal price (t) * Abnahmenahmemenge (t)) / Abnahmemenge.sum() for each bus
        # Abnahmemenge
        # LCOH, LCOE // weighted, not weighted // export, no export, mixed

        # Get hydrogen buses
        buses_h_export = n.buses.index[n.buses.index == "H2 export bus"]
        buses_h_noexport = n.buses.index[(n.buses.index.str[-2:] == "H2")]
        buses_h_mixed = buses_h_noexport.append(buses_h_export)

        # Get electricity buses (no differentiation between export and no export, because there is no export bus for electricity)
        buses_e = n.buses.index[(n.buses.index.str[-2:] == "AC")]

        # Calculate the not weighted ("notw") marginal price of the hydrogen buses
        lcoh_notw_export = (
            n.buses_t.marginal_price[buses_h_export].mean().mean().round(2)
        )
        lcoh_notw_noexport = (
            n.buses_t.marginal_price[buses_h_noexport].mean().mean().round(2)
        )
        lcoh_notw_mixed = n.buses_t.marginal_price[buses_h_mixed].mean().mean().round(2)

        # Calculate the not weighted ("notw") marginal price of the electricity buses
        lcoe_notw = n.buses_t.marginal_price[buses_e].mean().mean().round(2)

        # Calculate the weighted ("w") marginal price of the hydrogen buses

        # busname = "H2 export bus"
        # busname = ["MA.1.2_1_AC H2", "MA.1.3_1_AC H2"]
        # demand = get_bus_demand(n, busname)

        demand_h_export = get_bus_demand(n, buses_h_export)
        demand_h_noexport = get_bus_demand(n, buses_h_noexport)
        demand_h_mixed = get_bus_demand(n, buses_h_mixed)

        # Get the demand (load, links, stores?) of the hydrogen export bus

        load_h_export = n.loads_t.p["H2 export load"].sum()

        # PyPSA-Eur Approach
        # def calculate_prices(n, label, prices):
        #     prices = prices.reindex(prices.index.union(n.buses.carrier.unique()))

        #     # WARNING: this is time-averaged, see weighted_prices for load-weighted average
        #     prices[label] = n.buses_t.marginal_price.mean().groupby(n.buses.carrier).mean()

        #     return prices

        ############################################################

        # Save the cost and lcoh in the array according to the h2export and opts values using concat function
        cost_df = pd.concat(
            [
                cost_df,
                pd.DataFrame(
                    {
                        "h2export": [h2export],
                        "opts": [opts],
                        "cost": [cost],
                        "lcoh": [lcoh.values[0]],
                        "lcoh_marginal": [lcoh_marginal],
                        "lcoh_marginal_export": [lcoh_marginal_export],
                        "lcoe_marginal": [lcoe_marginal],
                        "lcoh_marginal_weighted": [lcoh_marginal_weighted],
                    }
                ),
            ],
            ignore_index=True,
        )

    print(cost0_df)
    print(cost_df)

    # Save the cost file
    cost_df.to_csv(snakemake.output.stats, index=False)

    # Plot the cost of the networks
    # cost_df.pivot(index='h2export', columns='opts', values='cost').plot(kind='bar')
    # plt.ylabel('Cost [EUR]')
    # plt.xlabel('H2 export [TWh]')
    # plt.title('Cost of the networks')
    # plt.savefig(snakemake.output.cost_plot)
