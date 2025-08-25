import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_visualization']['base_path']
subsets_filtered_dir_rel = config['dim_reduction_visualization']['subsets_filtered_dir']
umap_qt_dir_rel = config['dim_reduction_visualization']['umap_qt_dir']
umap_qt_coloring_br_tr_dir_rel = config['dim_reduction_visualization']['umap_qt_coloring_br_tr_dir']

# Construct full paths
subsamples_dir = os.path.join(base_path, subsets_filtered_dir_rel)
embeddings_dir = os.path.join(base_path, umap_qt_dir_rel)
output_directory = os.path.join(base_path, umap_qt_coloring_br_tr_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Define color palettes for Tech_replica and Metadata_Biorep
tech_replica_colors = {
    'B': '#2c7bb6',
    'C': '#fdae61',
    'D': '#d7191c'
}

biorep_colors = {
    'BR1': '#2c7bb6',
    'BR2': '#abd9e9',
    'BR3': '#fdae61',
    'BR4': '#d7191c'
}

# Function to create UMAP plots
def create_umap_plots(plot_data, color_dict, metadata_column, base_file_name):
    # Generate the full plot
    plt.figure(figsize=(10, 8))
    full_plot = sns.scatterplot(
        x='UMAP-1',
        y='UMAP-2',
        hue=metadata_column,
        data=plot_data,
        palette=color_dict,
        s=9,
        hue_order=sorted(color_dict.keys()),
        legend='brief'
    )
    plt.title(f'UMAP plot colored by {metadata_column}')
    plt.xlabel('UMAP-1')
    plt.ylabel('UMAP-2')
    plt.legend(title=metadata_column, loc='best', bbox_to_anchor=None)
    plt.savefig(os.path.join(output_directory, f"umap_{base_file_name}_{metadata_column}_full.png"), dpi=300)
    plt.close()

    # Create individual plots for each value
    for value in color_dict.keys():
        plot_data['Color'] = np.where(plot_data[metadata_column] == value, f'{metadata_column} {value}', 'Other')
        plot_data['Size'] = np.where(plot_data[metadata_column] == value, 20, 9)

        plt.figure(figsize=(10, 8))
        sns.scatterplot(
            x='UMAP-1',
            y='UMAP-2',
            data=plot_data[plot_data['Color'] == 'Other'],
            color='#E6E6E6',
            s=9,
            legend=False
        )
        sns.scatterplot(
            x='UMAP-1',
            y='UMAP-2',
            data=plot_data[plot_data['Color'] != 'Other'],
            color=color_dict[value],
            s=20,
            edgecolor='black',
            linewidth=0.2,
            legend=False
        )
        plt.title(f'UMAP plot colored by {metadata_column} {value}')
        plt.xlabel('UMAP-1')
        plt.ylabel('UMAP-2')

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[value], 
                   markersize=10, markeredgecolor='black', markeredgewidth=0.2, label=f'{metadata_column} {value}'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#E6E6E6', 
                   markersize=7, label='Other')
        ]
        plt.legend(handles=legend_elements, loc='best')
        #plt.gca().get_legend().remove()

        plt.savefig(os.path.join(output_directory, f"umap_{base_file_name}_{metadata_column}_{value}.png"), dpi=300)
        plt.close()

# List all the subsample files
subsample_files = [f for f in os.listdir(subsamples_dir) if f.endswith('.txt') and f != "feature_filtering_report.txt"]

# Loop through each subsample file
for subsample_file in subsample_files:
    file_path = os.path.join(subsamples_dir, subsample_file)
    data = pd.read_csv(file_path, sep='\t')
    
    # Randomize the order of rows in the dataset
    data = data.sample(frac=1, random_state=42).reset_index(drop=True)
    
    # Separate features and metadata
    metadata_columns = ['Concentration', 'counts_Cells', 'counts_Cytoplasm', 'counts_FilteredNuclei',
                        'Metadata_Well', 'Metadata_Day', 'Metadata_Biorep', 'Tech_replica', 'Day_Well_BR', 'cell_ID']
    features = data.drop(columns=metadata_columns)
    metadata = data[metadata_columns]
    
    # Replace infinite values with NaN and impute missing values
    features = features.replace([np.inf, -np.inf], np.nan)
    imputer = SimpleImputer(strategy='median')
    features_imputed = imputer.fit_transform(features)
    scaler = StandardScaler(with_mean=True, with_std=False)
    features_centered = scaler.fit_transform(features_imputed)
    
    # Load the saved UMAP embedding
    embedding_file_name = subsample_file.replace('.txt', '_umap_embedding_QuantileTransformer.npy')
    embedding_file_path = os.path.join(embeddings_dir, embedding_file_name)
    umap_embedding = np.load(embedding_file_path)
    
    # Prepare the plot data
    plot_data = pd.DataFrame(umap_embedding, columns=['UMAP-1', 'UMAP-2'])
    # Group E with B, F with C, and G with D
    plot_data['Tech_replica'] = metadata['Tech_replica'].replace({'E': 'B', 'F': 'C', 'G': 'D'})
    plot_data['Metadata_Biorep'] = metadata['Metadata_Biorep']
    
    # Determine the base file name without extension and directory
    base_file_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # Create plots for Tech_replica
    create_umap_plots(plot_data, tech_replica_colors, 'Tech_replica', base_file_name)
    
    # Create plots for Metadata_Biorep
    create_umap_plots(plot_data, biorep_colors, 'Metadata_Biorep', base_file_name)

print("UMAP plots using precomputed embeddings created for all subsamples!")