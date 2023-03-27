"""Script to compare networks.
"""

import pypsa
import pandas as pd
from matplotlib.pyplot import legend
import matplotlib.pyplot as plt
import numpy as np
import os
import logging


def plot_total_system_cost(n1, n2, n3):    

    n_names = ["Integrated", "Islanded"]
    costs_mar_es_only = [0, n1.objective/1e9]
    costs_mar_es_export = [n2.objective/1e9, 0]
    costs_export_only = [0, n3.objective/1e9]

    diff_percent = ((costs_export_only[1]+costs_mar_es_only[1])-costs_mar_es_export[0])*100/costs_mar_es_export[0]
    #bars = np.add(costs13, costs2).tolist()

    # The position of the bars on the x-axis
    r = [0,1]
    barWidth = 0.4

    # Create bars and stack them on previous bars
    plt.figure(figsize=(4, 3))
    plt.bar(r, costs_mar_es_only, color=['g'], alpha=1, edgecolor='white', width=barWidth)
    plt.bar(r, costs_mar_es_export, color=['g'], alpha=0.3, edgecolor='white', width=barWidth)
    plt.bar(r, costs_export_only, bottom=costs_mar_es_only, alpha=0.6, color=['g'], edgecolor='white', width=barWidth)

    # enhance graph
    plt.xticks(r, n_names)
    plt.ylabel("System cost in B€")
    plt.title(f"Saving {diff_percent:.1f} % of Total System Cost ({snakemake.wildcards.h2export} TWh export)")
    #plt.ylim(212, 214)

    # Add a legend
    plt.legend(["MAR ES only", "MAR ES + Export", "Export only"], loc='upper left', bbox_to_anchor=(1,1), ncol=1)

    # Add grid on y-axis
    plt.grid(axis='y', alpha=0.5)

    # Show graphic
    plt.savefig(snakemake.output.graph_systemcost, bbox_inches='tight')
    plt.show()

    # Print the costs costs_mar_es_only, costs_mar_es_export, costs_export_only
    print(f"Costs for MAR ES only: {costs_mar_es_only}")
    print(f"Costs for MAR ES + Export: {costs_mar_es_export}")
    print(f"Costs for Export only: {costs_export_only}")

    return




if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake('compare_networks', h2export=20)

        sets_path_to_root('aldehyde')

    n1 = pypsa.Network(snakemake.input.network_mar_es_only)
    n2 = pypsa.Network(snakemake.input.network_mar_es_export)
    n3 = pypsa.Network(snakemake.input.network_export_only)


    plot_total_system_cost(n1, n2, n3)
