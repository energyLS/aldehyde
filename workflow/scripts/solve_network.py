import pypsa
import os
import pandas as pd

def solve_network(n):


    n.lopf(pyomo=False, solver_name='gurobi')
    print(f"Solving complete with objective: {n.objective:.2f}")

    return



if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "solve_network",
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


    # Read network
    # https://pypsa.readthedocs.io/en/latest/components.html?highlight=override_component_attrs#custom-components
    n = pypsa.Network(snakemake.input.network) #, override_component_attrs=overrides)

    snapshots = n.snapshots

    solve_network(n)

    if pd.isna(n.objective):
        print("No solution found. Exiting.")
        exit()

    n.export_to_netcdf(snakemake.output[0])