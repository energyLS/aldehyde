"""Script to compare networks.
"""

import pypsa
import pandas as pd
from matplotlib.pyplot import legend
import matplotlib.pyplot as plt
import numpy as np
import os
import logging


def get_bus_demand(n, busname, carrier_limit=False, carrier_limit_integration=False):
    """Get the demand at a certain bus (includes stores/StorageUnits) based on the energy_balance() function. The demand excludes lines for electricity transport but includes H2 pipelines. It further excludes Stores and StorageUnits.

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

    # Filter for negative values, filter out Stores/StorageUnits, filter out export links "H2" which are not in "H2 export buses", and sum them up. Inclues H2 pipeline, excludes electricity lines. 
    demand = energy_balance_bus[energy_balance_bus < 0]
    demand = demand[(demand.index.get_level_values(0) != "Store") & (demand.index.get_level_values(0) != "StorageUnit")]
    demand = demand[(demand.index.get_level_values(0) != "Line")]
    #only when demand.index.get_level_values(1) == "H2 export bus", then H2 is excluded to avoid double counting through pipelines
    demand = demand[(demand.index.get_level_values(1) != "H2") | (demand.index.get_level_values(3) == "H2 export bus")]
    demand = demand.groupby("bus").sum()


    return demand.transpose()

def get_mean_prices(n):
    return (
        n.buses_t.marginal_price.mean()
        .groupby([n.buses.carrier])
        .first()
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("compare_integrated")

        sets_path_to_root("aldehyde")

    # Load networks from snakefile input
    # Create loop over snakemake inputs to load all networks

    # Loop over the networks and save the objective of the networks in array according to the h2export and opts wildcards
    metrics_df = pd.DataFrame(columns=["h2export", "opts", "cost"])
    cost0_df = pd.DataFrame(columns=["h2export", "opts", "cost"])

    for i in snakemake.input.networks:  # [0:limit]:
        n = pypsa.Network(i)
        cost = n.objective

        # Get the h2export (before "export") and opts values from the network name
        h2export = i.split("_")[-1][:-9]
        opts = i.split("_")[-6]
        statistics = n.statistics()
        energy_balance = n.statistics.energy_balance(aggregate_bus=False, aggregate_time=False)

        # Read and save the objective of the network if h2export is 0. Save the costs in a new array. Only works, if scenarios with h2export = 0 exist
        # if h2export == "0":
        #     n0 = pypsa.Network(i)
        #     cost0 = n0.objective
        #     # Save the cost0 in an array
        #     cost0_df = pd.concat(
        #         [
        #             cost0_df,
        #             pd.DataFrame(
        #                 {"h2export": [h2export], "opts": [opts], "cost": [cost0]}
        #             ),
        #         ],
        #         ignore_index=True,
        #     )

        # # Calculate the cost difference between the networks n and n0
        # cost_difference = cost - cost0_df[cost0_df["opts"] == opts]["cost"]  # in €

        # # Calculate the cost of the H2 export by taking the difference in system costs and dividing by the export demand
        # lcoh = cost_difference / (
        #     n.loads_t.p.loc[:, "H2 export load"].sum()
        #     * n.snapshot_weightings.generators[0]
        # )  # Cost difference between systems divided by the export demand


        ############################################################
        # New approach to calculate the weighted marginal price of the hydrogen buses: sum(Marginal price (t) * demand (t)) / demand.sum() for each bus. Then, the weighted mean (based on demand again) of all buses is calculated.
        # Step 1: get buses
        # Step 2: get demand at buses
        # Step 3: get (weighted) marginal price at buses

        
        carriers = ["H2", "AC", "oil", "gas"] # Get the buses with these carriers
        weighting_options = {"H2": [["Fischer-Tropsch", "inclusive", "all"],
                                    ["Fischer-Tropsch", "exclusive", "all"],
                                    [False, False, "exportonly"],
                                    [False, False, "noexport"],
                                    [False, False, "all"],
                                    ],
                             "AC": [["H2 Electrolysis", "exclusive", "all"], 
                                    ["H2 Electrolysis", "inclusive", "all"],
                                    [False, False, "all"]],
                             "oil": [[False, False, "all"]],
                             "gas": [[False, False, "all"]],
                            }
        marginals_df = pd.DataFrame(columns=[])
        expense_df = pd.DataFrame(columns=[])
        for carrier in carriers:

            for w_opts in weighting_options[carrier]:
                # Step 1: get buses with differentation for hydrogen export node
                if (carrier == "H2") & (w_opts[2] == "exportonly"):
                    buses_sel = n.buses.index[n.buses.index == "H2 export bus"]
                elif (carrier == "H2") & (w_opts[2] == "noexport"):
                    buses_sel = n.buses.index[(n.buses.index.str[-2:] == "H2")]
                else:
                    buses_sel = n.buses[n.buses.carrier == carrier].index

                # Step 2: get demand at buses
                demand = get_bus_demand(n, buses_sel, w_opts[0], w_opts[1])

                # Step 3: get (weighted) nodal marginal price at buses
                marginals_nodal = (n.buses_t.marginal_price[buses_sel] * demand).sum() / demand.sum()

                # Step 4.1: Nodal weighted mean of marginal price
                marginals = (marginals_nodal * demand.sum() / demand.sum().sum()).sum().round(2)

                # Step 4.2: Expense on carrier by consumer
                expense = (n.buses_t.marginal_price[buses_sel] * demand).sum().sum() * n.snapshot_weightings.generators[0] *(-1) /1e6 # in Mio. €

                # Get variable name  and save the marginal price and expense in the dataframe
                marginal_name = "mg_" + carrier + "_" + str(w_opts[1])[:5] + "_" + str(w_opts[0])[:5] + "_" + str(w_opts[2])
                expense_name = "exp_" + carrier + "_" + str(w_opts[1])[:5] + "_" + str(w_opts[0])[:5] + "_" + str(w_opts[2])
                marginals_df[marginal_name] = [marginals]
                expense_df[expense_name] = [expense]


##########################################
    
        # Get marginal prices of all buses, unweighted
        marginals_means = get_mean_prices(n)      
        marginals_means_df = pd.DataFrame(marginals_means).transpose()
        marginals_means_df = marginals_means_df.add_prefix("mean_mg_")


        # Get storage capacities
        H2_GWh = n.stores[(n.stores.carrier=="H2") & (n.stores.bus != "H2 export bus")].e_nom_opt.sum() / 1e3 # in GWh
        Battery_GWh = n.stores[(n.stores.carrier=="battery")].e_nom_opt.sum() / 1e3 # in GWh
        EV_Battery_GWh = n.stores[(n.stores.carrier=="Li ion")].e_nom_opt.sum() / 1e3 # in GWh
        H2export_GWh = n.stores[(n.stores.carrier=="H2") & (n.stores.bus == "H2 export bus")].e_nom_opt.sum() / 1e3 # in GWh
        ratio_H2_Battery = (H2_GWh + H2export_GWh) / Battery_GWh

        #H2_GWh, Battery_GWh, H2export_GWh = get_storage_capacities(n)

        # Get curtailment rates
        curtailmentrate_solar = statistics.loc["Generator", "Solar"].Curtailment / statistics.loc["Generator", "Solar"].Dispatch *100
        curtailmentrate_wind = statistics.loc["Generator", "Onshore Wind"].Curtailment / statistics.loc["Generator", "Onshore Wind"].Dispatch *100

        # Get base electricity demand
         # Get all loads which are connected to the bus with carrier "AC"
        loads = n.loads_t.p.groupby(n.loads.carrier, axis=1).sum()
        ac_buses = n.buses[n.buses.carrier == "AC"]
        ac_loads = n.loads.loc[n.loads.bus.isin(ac_buses.index)].carrier.unique()
        # Sum of demand through loads on AC buses
        el_base_demand = (loads.loc[:, ac_loads].sum() * n.snapshot_weightings.generators[0] / 1e6).sum().round(2) # TWh


        ### Calculate hydrogen cost composition (CAPEX, OPEX, electricity cost)
        # Get data from statistics
        capex = statistics.loc["Link", "H2 Electrolysis"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        opex = statistics.loc["Link", "H2 Electrolysis"].loc["Operational Expenditure"].round(2) / 1e6 # in Mio. €
        supply = statistics.loc["Link", "H2 Electrolysis"].loc["Supply"].sum() / 1e6 # in TWh

        # Determine expense on electricity for hydrogen production
        buses_sel = n.buses[n.buses.carrier == "AC"].index
        prices = n.buses_t.marginal_price[buses_sel] 
        demand = n.links_t.p1[n.links[n.links.carrier=="H2 Electrolysis"].index]
        demand.columns = prices.columns # Adjust column names to match prices for multiplication
        e_cost = (prices * demand).sum().sum() * (-1) * n.snapshot_weightings.generators[0] / 1e6 # in Mio. €

        lcoh_compo = ((capex + opex + e_cost) / supply).round(2) # in €/MWh
        capex_share = (capex / (capex + opex + e_cost) * 100).round(2)
        capex_ely = capex
        capex_ely_rel = capex_ely / supply
        opex_ely = opex + e_cost
        opex_ely_rel = opex_ely / supply

        # Marginal price of co2
        mg_co2 = n.buses_t.marginal_price.loc[:, n.buses[n.buses.carrier == "co2"].index].mean()[0].round(2)

        # Get capacities and CAPEX in generation technologies
        pv_capex = statistics.loc["Generator", "Solar"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        pv_p_nom_opt = statistics.loc["Generator", "Solar"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        pv_supply = statistics.loc["Generator", "Solar"].loc["Supply"].sum() / 1e6 # in TWh
        pv_cf = statistics.loc["Generator", "Solar"].loc["Capacity Factor"].round(2)

        onshore_capex = statistics.loc["Generator", "Onshore Wind"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        onshore_p_nom_opt = statistics.loc["Generator", "Onshore Wind"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        onshore_supply = statistics.loc["Generator", "Onshore Wind"].loc["Supply"].sum() / 1e6 # in TWh
        onshore_cf = statistics.loc["Generator", "Onshore Wind"].loc["Capacity Factor"].round(2)

        coal_capex = statistics.loc["Generator", "Coal"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        coal_p_nom_opt = statistics.loc["Generator", "Coal"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        coal_supply = statistics.loc["Generator", "Coal"].loc["Supply"].sum() / 1e6 # in TWh
        coal_cf = statistics.loc["Generator", "Coal"].loc["Capacity Factor"].round(2)

        ccgt_capex = statistics.loc["Generator", "Combined-Cycle Gas"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        ccgt_p_nom_opt = statistics.loc["Generator", "Combined-Cycle Gas"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        ccgt_supply = statistics.loc["Generator", "Combined-Cycle Gas"].loc["Supply"].sum() / 1e6 # in TWh
        ccgt_cf = statistics.loc["Generator", "Combined-Cycle Gas"].loc["Capacity Factor"].round(2)

        ror_capex = statistics.loc["Generator", "Run of River"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        ror_p_nom_opt = statistics.loc["Generator", "Run of River"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        ror_supply = statistics.loc["Generator", "Run of River"].loc["Supply"].sum() / 1e6 # in TWh
        ror_cf = statistics.loc["Generator", "Run of River"].loc["Capacity Factor"].round(2)

        oil_capex = statistics.loc["Generator", "Oil"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        oil_p_nom_opt = statistics.loc["Generator", "Oil"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        oil_supply = statistics.loc["Generator", "Oil"].loc["Supply"].sum() / 1e6 # in TWh
        oil_cf = statistics.loc["Generator", "Oil"].loc["Capacity Factor"].round(2)
        
        ocgt_capex = statistics.loc["Link", "Open-Cycle Gas"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        ocgt_p_nom_opt = statistics.loc["Link", "Open-Cycle Gas"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        ocgt_supply = statistics.loc["Link", "Open-Cycle Gas"].loc["Supply"].sum() / 1e6 # in TWh
        ocgt_cf = statistics.loc["Link", "Open-Cycle Gas"].loc["Capacity Factor"].round(2)

        ft_capex = statistics.loc["Link", "Fischer-Tropsch"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €
        ft_p_nom_opt = statistics.loc["Link", "Fischer-Tropsch"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        ft_supply = statistics.loc["Link", "Fischer-Tropsch"].loc["Supply"].sum() / 1e6 # in TWh
        ft_cf = statistics.loc["Link", "Fischer-Tropsch"].loc["Capacity Factor"].round(2)

        # Capacity factor electrolysis
        cf_electrolysis = statistics.loc["Link", "H2 Electrolysis"].loc["Capacity Factor"].round(2)
        electrolysis_p_nom_opt = statistics.loc["Link", "H2 Electrolysis"].loc["Optimal Capacity"].round(2) / 1e3 # in GW
        electrolysis_supply = statistics.loc["Link", "H2 Electrolysis"].loc["Supply"].sum() / 1e6 # in TWh
        electrolysis_capex = statistics.loc["Link", "H2 Electrolysis"].loc["Capital Expenditure"].round(2) / 1e6 # in Mio. €


        # Save the cost and lcoh in the array according to the h2export and opts values using concat function
        metrics_df = pd.concat(
            [
                metrics_df,
                pd.concat([
                    marginals_df,
                    expense_df,
                    marginals_means_df,
                    pd.DataFrame(
                        {
                            "h2export": [h2export],
                            "opts": [opts],
                            "cost": [cost],
                            "mg_co2": [mg_co2],
                            #"lcoh_system": [lcoh.values[0]],
                            "lcoh_compo": [lcoh_compo],
                            "capex_share": [capex_share],
                            "capex_ely": [capex_ely],
                            "capex_ely_rel": [capex_ely_rel],
                            "opex_ely": [opex_ely],
                            "opex_ely_rel": [opex_ely_rel],
                            "H2_GWh": [H2_GWh],
                            "Battery_GWh": [Battery_GWh],
                            "EV_Battery_GWh": [EV_Battery_GWh],
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
                            "pv_supply": [pv_supply],
                            "pv_cf": [pv_cf],
                            "onshore_supply": [onshore_supply],
                            "onshore_cf": [onshore_cf],
                            "coal_supply": [coal_supply],
                            "coal_cf": [coal_cf],
                            "ccgt_supply": [ccgt_supply],
                            "ccgt_cf": [ccgt_cf],
                            "ror_supply": [ror_supply],
                            "ror_cf": [ror_cf],
                            "oil_supply": [oil_supply],
                            "oil_cf": [oil_cf],
                            "ocgt_supply": [ocgt_supply],
                            "ocgt_cf": [ocgt_cf],
                            "ft_capex": [ft_capex],
                            "ft_p_nom_opt": [ft_p_nom_opt],
                            "ft_supply": [ft_supply],
                            "ft_cf": [ft_cf],
                            "electrolysis_p_nom_opt": [electrolysis_p_nom_opt],
                            "electrolysis_supply": [electrolysis_supply],
                            "electrolysis_capex": [electrolysis_capex],
                        }
                    ),
                ], axis=1),
            ],
            ignore_index=True,
        )

    # Save the cost file
    metrics_df.to_csv(snakemake.output.stats, index=False)