"""Script to compare networks.
"""

import pypsa
import pandas as pd
from matplotlib.pyplot import legend
import matplotlib.pyplot as plt
import numpy as np
import os
import logging


def get_buses(n):
    # Get hydrogen buses
    buses_h_export = n.buses.index[n.buses.index == "H2 export bus"]
    buses_h_noexport = n.buses.index[(n.buses.index.str[-2:] == "H2")]
    buses_h_mixed = buses_h_noexport.append(buses_h_export)

    # Get electricity buses (no differentiation between export and no export, because there is no export bus for electricity)
    buses_e = n.buses.index[(n.buses.index.str[-2:] == "AC")]
    return buses_h_export, buses_h_noexport, buses_h_mixed, buses_e


def get_bus_demand(n, busname, carrier_limit=False, carrier_limit_integration=False):
    """Get the demand at a certain bus (includes stores/StorageUnits) based on the energy_balance() function

    Parameters
    ----------
    busname : _type_
        _description_
    n : _type_
        _description_
    """

    # Get the energy balance of the bus specified in busname
    if carrier_limit == False:
        energy_balance_bus = energy_balance.loc[:, :, :, busname]
    else:
        if carrier_limit_integration == "inclusive":
            energy_balance_bus = energy_balance.loc[:,carrier_limit, :, busname]

        elif carrier_limit_integration == "exclusive":  # exclude carrier_limit
            energy_balance_bus = energy_balance.loc[:, :, :, busname]
            # Create a boolean mask to exclude rows with 'carrier_limit' in the third level
            condition = energy_balance_bus.index.get_level_values(1) != carrier_limit
            energy_balance_bus = energy_balance_bus[condition]  
        else:
            logging.error("carrier_limit_integration must be 'inclusive' or 'exclusive'")

    # Filter for negative values and sum them up (note: may include stores/StorageUnits)
    demand = energy_balance_bus[energy_balance_bus < 0]
    demand = demand.groupby("bus").sum()

    return demand.transpose()


def calc_notw_marginals(n, buses_h_export, buses_h_noexport, buses_h_mixed, buses_e):
    # Calculate the not weighted ("notw") marginal price of the hydrogen buses
    lcoh_notw_export = n.buses_t.marginal_price[buses_h_export].mean().mean().round(2)
    lcoh_notw_noexport = (
        n.buses_t.marginal_price[buses_h_noexport].mean().mean().round(2)
    )
    lcoh_notw_mixed = n.buses_t.marginal_price[buses_h_mixed].mean().mean().round(2)

    # Calculate the not weighted ("notw") marginal price of the electricity buses
    lcoe_notw = n.buses_t.marginal_price[buses_e].mean().mean().round(2)

    return lcoh_notw_export, lcoh_notw_noexport, lcoh_notw_mixed, lcoe_notw


def calc_weighted_marginals(
    n, buses_h_export, buses_h_noexport, buses_h_mixed, buses_e
):
    # Get the nodal weighted marginal price of the hydrogen buses: sum(Marginal price (t) * Abnahmenahmemenge (t)) / Abnahmemenge.sum() for each bus
    lcoh_w_noexport_nodal = (
        n.buses_t.marginal_price[buses_h_noexport] * demand_h_noexport
    ).sum() / demand_h_noexport.sum()
    lcoh_w_export_nodal = (
        n.buses_t.marginal_price[buses_h_export] * demand_h_export
    ).sum() / demand_h_export.sum()
    lcoh_w_mixed_nodal = (
        n.buses_t.marginal_price[buses_h_mixed] * demand_h_mixed
    ).sum() / demand_h_mixed.sum()

    # Also for electricity buses
    lcoe_w_nodal = (n.buses_t.marginal_price[buses_e] * demand_e).sum() / demand_e.sum()
    lcoe_w_electrolysis_nodal = (
        n.buses_t.marginal_price[buses_e] * demand_e_electrolysis
    ).sum() / demand_e_electrolysis.sum()
    lcoe_w_no_electrolysis_nodal = (
        n.buses_t.marginal_price[buses_e] * demand_e_no_electrolysis
    ).sum() / demand_e_no_electrolysis.sum()

    # Get the weighted mean of the nodal weighted marginal price of the hydrogen buses
    lcoh_w_noexport = (
        (
            lcoh_w_noexport_nodal
            * demand_h_noexport.sum()
            / demand_h_noexport.sum().sum()
        )
        .sum()
        .round(2)
    )
    lcoh_w_export = (
        (lcoh_w_export_nodal * demand_h_export.sum() / demand_h_export.sum().sum())
        .sum()
        .round(2)
    )
    lcoh_w_mixed = (
        (lcoh_w_mixed_nodal * demand_h_mixed.sum() / demand_h_mixed.sum().sum())
        .sum()
        .round(2)
    )

    # Also for electricity buses
    lcoe_w = (lcoe_w_nodal * demand_e.sum() / demand_e.sum().sum()).sum().round(2)
    lcoe_w_electrolysis = (
        lcoe_w_electrolysis_nodal * demand_e_electrolysis.sum() / demand_e_electrolysis.sum().sum()
    ).sum().round(2)
    lcoe_w_no_electrolysis = (
        lcoe_w_no_electrolysis_nodal * demand_e_no_electrolysis.sum() / demand_e_no_electrolysis.sum().sum()
    ).sum().round(2)

    return lcoh_w_noexport, lcoh_w_export, lcoh_w_mixed, lcoe_w, lcoe_w_electrolysis, lcoe_w_no_electrolysis


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
    metrics_df = pd.DataFrame(columns=["h2export", "opts", "cost", "lcoh"])
    cost0_df = pd.DataFrame(columns=["h2export", "opts", "cost"])

    for i in snakemake.input.networks:  # [0:limit]:
        n = pypsa.Network(i)
        cost = n.objective

        # Get the h2export (before "export") and opts values from the network name
        h2export = i.split("_")[-1][:-9]
        opts = i.split("_")[-6]
        statistics = n.statistics()
        energy_balance = n.statistics.energy_balance(aggregate_bus=False, aggregate_time=False)

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


        ############################################################
        # New approach to calculate the weighted marginal price of the hydrogen buses: sum(Marginal price (t) * demand (t)) / demand.sum() for each bus. Then, the weighted mean (based on demand again) of all buses is calculated.
        # Abnahmemenge
        # LCOH, LCOE // weighted, not weighted // export, no export, mixed

        buses_h_export, buses_h_noexport, buses_h_mixed, buses_e = get_buses(n)

        (
            lcoh_notw_export,
            lcoh_notw_noexport,
            lcoh_notw_mixed,
            lcoe_notw,
        ) = calc_notw_marginals(
            n, buses_h_export, buses_h_noexport, buses_h_mixed, buses_e
        )

        # Calculate the weighted ("w") marginal price of the hydrogen buses
        demand_h_export = get_bus_demand(n, buses_h_export)
        demand_h_noexport = get_bus_demand(n, buses_h_noexport)
        demand_h_mixed = get_bus_demand(n, buses_h_mixed)
        demand_e = get_bus_demand(n, buses_e)

        demand_e_electrolysis = get_bus_demand(n, buses_e, "H2 Electrolysis", "inclusive")
        demand_e_no_electrolysis = get_bus_demand(n, buses_e, "H2 Electrolysis", "exclusive")

        #lcoe_w = calc_weighted_marginals_new(

        lcoh_w_noexport, lcoh_w_export, lcoh_w_mixed, lcoe_w, lcoe_w_electrolysis, lcoe_w_no_electrolysis = calc_weighted_marginals(
            n, buses_h_export, buses_h_noexport, buses_h_mixed, buses_e
        )

        # Get storage capacities
        H2_GWh = n.stores[(n.stores.carrier=="H2") & (n.stores.bus != "H2 export bus")].e_nom_opt.sum() / 1e3 # in GWh
        Battery_GWh = n.stores[(n.stores.carrier=="battery")].e_nom_opt.sum() / 1e3 # in GWh
        H2export_GWh = n.stores[(n.stores.carrier=="H2") & (n.stores.bus == "H2 export bus")].e_nom_opt.sum() / 1e3 # in GWh
        ratio_H2_Battery = (H2_GWh + H2export_GWh) / Battery_GWh

        #H2_GWh, Battery_GWh, H2export_GWh = get_storage_capacities(n)

        # Get curtailment rates
        curtailmentrate_solar = statistics.loc["Generator", "Solar"].Curtailment / statistics.loc["Generator", "Solar"].Dispatch *100
        curtailmentrate_wind = statistics.loc["Generator", "Onshore Wind"].Curtailment / statistics.loc["Generator", "Onshore Wind"].Dispatch *100

        # Get base electricity demand
        el_base_demand = n.loads_t.p_set[n.loads[n.loads.carrier=="AC"].index].sum().sum()/1e6 * n.snapshot_weightings.generators[0]# in TWh

        # Capacity factor electrolysis
        cf_electrolysis = statistics.loc["Link", "H2 Electrolysis"].loc["Capacity Factor"].round(2)

        # Get capacities and CAPEX in generation technologies
        pv_capex = statistics.loc["Generator", "Solar"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        pv_p_nom_opt = statistics.loc["Generator", "Solar"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        onshore_capex = statistics.loc["Generator", "Onshore Wind"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        onshore_p_nom_opt = statistics.loc["Generator", "Onshore Wind"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        coal_capex = statistics.loc["Generator", "Coal"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        coal_p_nom_opt = statistics.loc["Generator", "Coal"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        ccgt_capex = statistics.loc["Generator", "Combined-Cycle Gas"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        ccgt_p_nom_opt = statistics.loc["Generator", "Combined-Cycle Gas"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        ror_capex = statistics.loc["Generator", "Run of River"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        ror_p_nom_opt = statistics.loc["Generator", "Run of River"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        oil_capex = statistics.loc["Generator", "Oil"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        oil_p_nom_opt = statistics.loc["Generator", "Oil"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        
        ocgt_capex = statistics.loc["Link", "Open-Cycle Gas"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        ocgt_p_nom_opt = statistics.loc["Link", "Open-Cycle Gas"].loc["Optimal Capacity"].round(2) / 1e3 # in GW


        # Save the cost and lcoh in the array according to the h2export and opts values using concat function
        metrics_df = pd.concat(
            [
                metrics_df,
                pd.DataFrame(
                    {
                        "h2export": [h2export],
                        "opts": [opts],
                        "cost": [cost],
                        "lcoh_system": [lcoh.values[0]],
                        "lcoh_notw_export": [lcoh_notw_export],
                        "lcoh_notw_noexport": [lcoh_notw_noexport],
                        "lcoh_notw_mixed": [lcoh_notw_mixed],
                        "lcoe_notw": [lcoe_notw],
                        "lcoh_w_export": [lcoh_w_export],
                        "lcoh_w_noexport": [lcoh_w_noexport],
                        "lcoh_w_mixed": [lcoh_w_mixed],
                        "lcoe_w": [lcoe_w],
                        "lcoe_w_electrolysis": [lcoe_w_electrolysis],
                        "lcoe_w_no_electrolysis": [lcoe_w_no_electrolysis],
                        "H2_GWh": [H2_GWh],
                        "Battery_GWh": [Battery_GWh],
                        "H2export_GWh": [H2export_GWh],
                        "ratio_H2_Battery": [ratio_H2_Battery],
                        "curtailmentrate_solar": [curtailmentrate_solar],
                        "curtailmentrate_wind": [curtailmentrate_wind],
                        "el_base_demand": [el_base_demand],
                        "cf_electrolysis": [cf_electrolysis],
                        "pv_capex": [pv_capex],
                        "pv_p_nom_opt": [pv_p_nom_opt],
                        "onshore_capex": [onshore_capex],
                        "onshore_p_nom_opt": [onshore_p_nom_opt],
                        "coal_capex": [coal_capex],
                        "coal_p_nom_opt": [coal_p_nom_opt],
                        "ccgt_capex": [ccgt_capex],
                        "ccgt_p_nom_opt": [ccgt_p_nom_opt],
                        "ror_capex": [ror_capex],
                        "ror_p_nom_opt": [ror_p_nom_opt],
                        "oil_capex": [oil_capex],
                        "oil_p_nom_opt": [oil_p_nom_opt],
                        "ocgt_capex": [ocgt_capex],
                        "ocgt_p_nom_opt": [ocgt_p_nom_opt],
                    }
                ),
            ],
            ignore_index=True,
        )


    # Save the cost file
    metrics_df.to_csv(snakemake.output.stats, index=False)