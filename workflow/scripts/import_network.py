import pypsa
import os







if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake('import_network')

        sets_path_to_root('aldehyde')


    # Read European and Moroccan pypsa networks
    # https://pypsa.readthedocs.io/en/latest/components.html?highlight=override_component_attrs#custom-components
    n_eur = pypsa.Network(snakemake.input.network) #, override_component_attrs=overrides)