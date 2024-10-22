# ALDEHYDE - locAL DEcarbonization and HYDrogen Export

This repository contains the entire scientific project, including code and report for the paper **"The impact of temporal hydrogen regulation on hydrogen exporters and their domestic energy transition"**.


## Abstract
As global demand for green hydrogen rises, potential hydrogen exporters move into the spotlight. However, the large-scale installation of on-grid hydrogen electrolysis for export can have profound impacts on domestic energy prices and energy-related emissions. Our investigation explores the interplay of hydrogen exports, domestic energy transition and temporal hydrogen regulation, employing a sector-coupled energy model in Morocco. We find substantial co-benefits of domestic climate change mitigation and hydrogen exports, whereby exports can reduce domestic electricity prices while mitigation reduces hydrogen export prices. However, increasing hydrogen exports quickly in a system that is still dominated by fossil fuels can substantially raise domestic electricity prices, if green hydrogen production is not regulated. Surprisingly, temporal matching of hydrogen production lowers domestic electricity cost by up to 31% while the effect on exporters is minimal. This policy instrument can steer the welfare (re-)distribution between hydrogen exporting firms, hydrogen importers, and domestic electricity consumers and hereby increases acceptance among actors.


## Installation and Usage

1. Open your terminal at a location where you want to install the repository aldehyde including it's subworkflows PyPSA-Earth and PyPSA-Earth-Sec. Type the following in your terminal to download the package and the dependency (pypsa-earth) from GitHub.
   Note that the tag `--recursive-submodules` is needed to automatically clone also the pypsa-earth dependency.

   ```bash
       .../some/path/without/spaces % git clone --recurse-submodules https://github.com/energyLS/aldehyde.git
   ```


2. Move the current directory to the head of the repository.
   ```bash
       .../some/path/without/spaces % cd aldehyde
   ```



4. The python package requirements are curated in the `workflow/subworkflows/pypsa-earth-sec/pypsa-earth/envs/environment.yaml` file of the pypsa-earth repository. The environment can be installed using `conda` or `mamba`:

   ```bash
        cd aldehyde/workflow/subworkflows/pypsa-earth-sec
       .../aldehyde/pypsa-earth-sec % conda env create -f pypsa-earth/envs/environment.yaml
   ```

5. For running the optimization one has to install the solver. We can recommend the open source HiGHs solver, see more details on solvers in the documentation of [PyPSA-Earth](https://pypsa-earth.readthedocs.io/en/latest/).

The total installation time of cloning the repository and installing the environment is approximately 30 mins, given the prior installation of conda or mamba.



## Repository structure

* `config`: contains configuration files for aldehyde (`config.yaml`) and PyPSA-Earth-Sec (`config.pypsa-earth-sec.yaml`) for high-level plotting

* `report`: contains the .tex files for the paper
* `workflow/notebooks`: contains the Jupyter notebooks used for the evaluation of results
* `workflow/scripts`: contains the scripts used for the evaluation of results
* `workflow/subworkflows`: contains the `PyPSA-Earth-Sec` workflow which includes the `PyPSA-Earth` workflow. PyPSA-Earth-Sec is based on the configuration in `config.paper.yaml` and PyPSA-Earth is based on the configuration in `config.pypsa-earth.yaml`.

## Run scenarios
For running the model, navigate to the PyPSA-Earth-Sec model by:
```bash
cd workflow/subworkflows/pypsa-earth-sec
```
To solve all networks, run the following command:
```bash
snakemake -j 1 solve_all_networks -n
```

Please follow the documentation of [PyPSA-Earth](https://pypsa-earth.readthedocs.io/en/latest/) and the [Readme of PyPSA-Earth-Sec](https://github.com/pypsa-meets-earth/pypsa-earth-sec/blob/main/README.md) for more details. The estimated time to run one single optimization is 40 mins on a standard laptop, the full set of paper results includes over 360 optimizations. To run the full set, a high-performance computer is recommended.

## Reproducibility
The paper results and analysis are created on the following commits:s

* `aldehyde`: on commit https://github.com/energyLS/aldehyde/commit/465750d6f12716c44f77980e8ea56f05997c20ba which includes the submodule of

* `PyPSA-Earth-Sec` on the commit https://github.com/pypsa-meets-earth/pypsa-earth-sec/tree/6ab3255d5b6f5f9182ddddc04da658ab1902f975 which includes the submodule of

* `PyPSA-Earth` on the commit https://github.com/pypsa-meets-earth/pypsa-earth/tree/84a0aa4470be9663657aa17540cdf08c8fa0f0b6

## Result data
A dataset of the model results is available on \href{https://doi.org/10.5281/zenodo.10951650}{Zenodo} under a CC-BY-4.0 license.

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
