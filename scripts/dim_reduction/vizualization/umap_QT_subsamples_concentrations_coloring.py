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
umap_qt_coloring_dir_rel = config['dim_reduction_visualization']['umap_qt_coloring_dir']

# Construct full paths
subsamples_dir = os.path.join(base_path, subsets_filtered_dir_rel)
embeddings_dir = os.path.join(base_path, umap_qt_dir_rel)
output_directory = os.path.join(base_path, umap_qt_coloring_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Define the color palette
colors = [
    "#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8",
    "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"
]

# Function to assign colors to concentrations
def assign_colors(concentrations):
    color_dict = {conc: colors[i % len(colors)] for i, conc in enumerate(sorted(concentrations, reverse=True))}
    return color_dict

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
    plot_data['Concentration'] = metadata['Concentration'].values
    
    # Create a color dictionary
    concentrations = sorted(plot_data['Concentration'].unique())
    color_dict = assign_colors(concentrations)
    
    # Generate the full plot and capture the color mapping
    plt.figure(figsize=(10, 8))
    full_plot = sns.scatterplot(
        x='UMAP-1',
        y='UMAP-2',
        hue='Concentration',
        data=plot_data,
        palette=color_dict,
        s=9,
        hue_order=sorted(concentrations, reverse=True),
        #legend='brief'
    )
    plt.title('UMAP plot colored by Concentration')
    plt.xlabel('UMAP-1')
    plt.ylabel('UMAP-2')
    
    # Adjust legend
    #plt.legend(title='Concentration', loc='best', bbox_to_anchor=None)
    plt.gca().get_legend().remove()
    
    # Determine the base file name without extension and directory
    base_file_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # Save the full plot
    plt.savefig(os.path.join(output_directory, f"umap_{base_file_name}_full_tc_fixed.png"), dpi=300)
    plt.close()
    
    # Loop through the specified concentrations and create a plot for each
    for concentration in concentrations:
        plot_data['Color'] = np.where(plot_data['Concentration'] == concentration, f'Concentration {concentration}', 'Other')
        plot_data['Size'] = np.where(plot_data['Concentration'] == concentration, 20, 9)
        
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
        
        # Plot highlighted concentration points with edge
        sns.scatterplot(
            x='UMAP-1',
            y='UMAP-2',
            data=plot_data[plot_data['Color'] != 'Other'],
            color=color_dict[concentration],
            s=20,
            edgecolor='black',
            linewidth=0.2,
            legend=False
        )
        
        plt.title(f'UMAP plot colored by Concentration {concentration}')
        plt.xlabel('UMAP-1')
        plt.ylabel('UMAP-2')
        
        # Create a custom legend
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[concentration], 
                   markersize=10, markeredgecolor='black', markeredgewidth=0.2, label=f'Concentration {concentration}'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#E6E6E6', 
                   markersize=7, label='Other')
        ]
        plt.legend(handles=legend_elements, loc='best')
    
        # Construct the output file paths
        png_file_path = os.path.join(output_directory, f"umap_{base_file_name}_c{concentration}.png")
    
        # Save the figures
        plt.savefig(png_file_path, dpi=300)
        plt.close()

print("UMAP plots using precomputed embeddings created for all subsamples!")