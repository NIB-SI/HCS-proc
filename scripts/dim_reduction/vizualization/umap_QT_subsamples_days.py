import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import QuantileTransformer
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_visualization']['base_path']
subsets_days_per_conc_filtered_dir_rel = config['dim_reduction_visualization']['subsets_days_per_conc_filtered_dir']
umap_qt_days_per_conc_dir_rel = config['dim_reduction_visualization']['umap_qt_days_per_conc_dir']

# Construct full paths
subsamples_dir = os.path.join(base_path, subsets_days_per_conc_filtered_dir_rel)
output_directory = os.path.join(base_path, umap_qt_days_per_conc_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Set a global random seed
np.random.seed(42)

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
    
    # Use a fixed random state for sampling
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
    scaler = QuantileTransformer(random_state=42)
    features_scaled = scaler.fit_transform(features_imputed)
    
    # Run UMAP with fixed random state
    umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
    umap_embedding = umap_model.fit_transform(features_scaled)
    
    # Prepare the plot data
    plot_data = pd.DataFrame(umap_embedding, columns=['UMAP-1', 'UMAP-2'])
    plot_data['Day'] = metadata['Metadata_Day'].values
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
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
    
    # Save the UMAP embedding for later use
    embedding_save_path = os.path.join(output_directory, f"{base_file_name}_umap_embedding_QuantileTransformer.npy")
    np.save(embedding_save_path, umap_embedding)
    
    # Construct the output file paths
    png_file_path = os.path.join(output_directory, f"umap_{base_file_name}_day_QuantileTransformer.png")
    
    # Save the figures
    plt.savefig(png_file_path, dpi=300)
    plt.close()

print("UMAP plots and embeddings created for all subsamples with QuantileTransformer scaling, colored by day!")