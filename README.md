# ALDEHYDE - locAL DEcarbonization and HYDrogen Export

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Run new scenarios
In `aldehyde`:

* in `config`: set `scenario`, `export`
    
* in `config.pypsa-earth-sec.yaml`: set `H2_network: True`, `H2_network_limit: 10000`, `run`

In `pypsa-earth-sec`:
* in `config.pypsa-earth.yaml`: set `["electricity"]["co2base"]: 40.0e+6`. Note: if you change that, make sure to delete the old networks
* in `config.pypsa-earth.yaml`: set `["cluster_options"]["alternative_clustering"]: True`

## Analyse new scenarios
In `aldehyde`:

* Contour -> `config.yaml`: takes the `scenario` and `export` params. Plot specifics in `plots["contour_plot"]`
* Spatials (e.g. potentials) -> `config.yaml`: takes the `plots["spatial_plot"]`
* PyPSA-Earth-Sec (e.g. Balances, H2-network) -> `config.pypsa-earth-sec.yaml`: takes the `scenario` and `export` params.



## Getting ready

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml

## Run the analysis

    snakemake -call

This will run all analysis steps to reproduce results and eventually build the report.

To generate a PDF of the dependency graph of all steps `build/dag.pdf` run:

    snakemake -c1 --use-conda -f dag


## Clone the repository

    git clone --recursive https://github.com/energyLS/aldehyde.git

## Config

Set the configuration in config/config.yaml. The scenario wildcards are propagated to all rules, so you can run the analysis for different scenarios.
Set PyPSA-Earth-Sec specific parameters in config/config.pypsa-eur-sec.yaml.

## Repo structure

* `config`: configurations used in the study
* `data`: place for raw data
* `report`: contains all files necessary to build the report; plots and result files are generated automatically
* `workflow`: contains the Snakemake workflow
* `build`: will contain all results (does not exist initially)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
