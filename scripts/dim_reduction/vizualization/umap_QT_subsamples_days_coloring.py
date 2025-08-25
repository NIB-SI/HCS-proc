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
subsets_days_per_conc_filtered_dir_rel = config['dim_reduction_visualization']['subsets_days_per_conc_filtered_dir']
umap_qt_days_per_conc_dir_rel = config['dim_reduction_visualization']['umap_qt_days_per_conc_dir']
umap_qt_days_coloring_dir_rel = config['dim_reduction_visualization']['umap_qt_days_coloring_dir']

# Construct full paths
subsamples_dir = os.path.join(base_path, subsets_days_per_conc_filtered_dir_rel)
embeddings_dir = os.path.join(base_path, umap_qt_days_per_conc_dir_rel)
output_directory = os.path.join(base_path, umap_qt_days_coloring_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Define the color palette for days
day_colors = {
    'D1': '#2c7bb6',
    'D5': '#abd9e9',
    'D7': '#fdae61',
    'D9': '#d7191c'
}

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
    plot_data['Day'] = metadata['Metadata_Day'].values
    
    # Generate the full plot and capture the color mapping
    plt.figure(figsize=(10, 8))
    full_plot = sns.scatterplot(
        x='UMAP-1',
        y='UMAP-2',
        hue='Day',
        data=plot_data,
        palette=day_colors,
        s=9,
        hue_order=['D1', 'D5', 'D7', 'D9'],
        #legend='brief'
    )
    plt.title('UMAP plot colored by Day')
    plt.xlabel('UMAP-1')
    plt.ylabel('UMAP-2')
    
    # Adjust legend
    #plt.legend(title='Day', loc='best', bbox_to_anchor=None)
    plt.gca().get_legend().remove()
    
    # Determine the base file name without extension and directory
    base_file_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # Save the full plot
    plt.savefig(os.path.join(output_directory, f"umap_{base_file_name}_full_tc_fixed.png"), dpi=300)
    plt.close()
    
    # Loop through the specified days and create a plot for each
    for day in ['D1', 'D5', 'D7', 'D9']:
        plot_data['Color'] = np.where(plot_data['Day'] == day, f'Day {day}', 'Other')
        plot_data['Size'] = np.where(plot_data['Day'] == day, 20, 9)
        
        # Create the plot
        plt.figure(figsize=(10, 8))
        
        # Plot background points first
        sns.scatterplot(
            x='UMAP-1',
            y='UMAP-2',
            data=plot_data[plot_data['Color'] == 'Other'],
            color='#E6E6E6',
            s=9,
            legend=False
        )
        
        # Plot highlighted day points with edge
        sns.scatterplot(
            x='UMAP-1',
            y='UMAP-2',
            data=plot_data[plot_data['Color'] != 'Other'],
            color=day_colors[day],
            s=20,
            edgecolor='black',
            linewidth=0.2,
            legend=False
        )
        
        plt.title(f'UMAP plot colored by Day {day}')
        plt.xlabel('UMAP-1')
        plt.ylabel('UMAP-2')
        
        # Create a custom legend
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=day_colors[day], 
                   markersize=10, markeredgecolor='black', markeredgewidth=0.2, label=f'Day {day}'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#E6E6E6', 
                   markersize=7, label='Other')
        ]
        plt.legend(handles=legend_elements, loc='best')
    
        # Construct the output file paths
        png_file_path = os.path.join(output_directory, f"umap_{base_file_name}_d{day}.png")
    
        # Save the figures
        plt.savefig(png_file_path, dpi=300)
        plt.close()

print("UMAP plots using precomputed embeddings created for all subsamples, colored by day!")