"""Script to compare networks.
"""

import pypsa
import pandas as pd
from matplotlib.pyplot import legend
import matplotlib.pyplot as plt
import numpy as np
import os
import logging






if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake('compare_integrated')

        sets_path_to_root('aldehyde')

    # Load networks from snakefile input
    # Create loop over snakemake inputs to load all networks

    # limit the number of networks
    #limit = 5
    

    # Loop over the networks and save the objective of the networks in array according to the h2export and sopts wildcards
    cost_df = pd.DataFrame(columns=['h2export', 'sopts', 'cost'])

    for i in snakemake.input: #[0:limit]:
        n = pypsa.Network(i)
        cost = n.objective
        
        # Get the h2export (before "export") and sopts values from the network name
        h2export = i.split('_')[-1][:-9]
        sopts = i.split('_')[-5]

        # Save the cost in the array according to the h2export and sopts values using concat function
        cost_df = pd.concat([cost_df, pd.DataFrame({'h2export': [h2export], 'sopts': [sopts], 'cost': [cost]})], ignore_index=True)

    print(cost_df)

    # Save the cost file
    cost_df.to_csv(snakemake.output.cost_file, index=False)

    # Plot the cost of the networks
    cost_df.pivot(index='h2export', columns='sopts', values='cost').plot(kind='bar')
    plt.ylabel('Cost [EUR]')
    plt.xlabel('H2 export [TWh]')
    plt.title('Cost of the networks')
    plt.savefig(snakemake.output.cost_plot)
    





