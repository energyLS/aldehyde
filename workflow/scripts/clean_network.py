import pypsa
import os
import pandas as pd
from _helpers import override_component_attrs

def remove_non_electric_buses(n):
    """
    remove buses from pypsa-eur with carriers which are not AC buses
    """
    print("drop buses from network with carrier: ", n.buses[~n.buses.carrier.isin(["AC", "DC"])].carrier.unique())
    n.buses = n.buses[n.buses.carrier.isin(["AC", "DC"])]

    return


def remove_elec_base_techs(n):
    """remove conventional generators (e.g. OCGT) and storage units (e.g. batteries and H2)
    from base electricity-only network, since they're added here differently using links
    """

    for c in n.iterate_components(snakemake.config["pypsa_earth"]):
        to_keep = snakemake.config["pypsa_earth"][c.name]
        to_remove = pd.Index(c.df.carrier.unique()).symmetric_difference(to_keep)
        print("Removing", c.list_name, "with carrier", to_remove)
        names = c.df.index[c.df.carrier.isin(to_remove)]
        n.mremove(c.name, names)
        n.carriers.drop(to_remove, inplace=True, errors="ignore")

# def manipulate_battery_capex(n, new_capex):

#     n.stores.loc[n.stores.carrier == "battery", "capital_cost"] = new_capex

#     return

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "clean_network",
            simpl="",
            clusters="10",
            ll="c1.0",
            opts="Co2L0.10",
            planning_horizons="2030",
            sopts="6H",
            discountrate=0.071,
            demand="DF",
            h2export=120,
        )

        sets_path_to_root('aldehyde')

    # Read network
    # https://pypsa.readthedocs.io/en/latest/components.html?highlight=override_component_attrs#custom-components
    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    remove_elec_base_techs(n)
    #remove_non_electric_buses(n)

    # Manipulate battery capex
    #new_capex = 10.0
    #manipulate_battery_capex(n, new_capex)

    # empty dataframe
    #n.global_constraints = n.global_constraints.iloc[0:0]

    

    n.export_to_netcdf(snakemake.output[0])