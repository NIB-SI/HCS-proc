import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import QuantileTransformer
from mpl_toolkits.mplot3d import Axes3D
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_visualization']['base_path']
subsets_filtered_dir_rel = config['dim_reduction_visualization']['subsets_filtered_dir']
umap_qt_3d_dir_rel = config['dim_reduction_visualization']['umap_qt_3d_dir']

# Construct full paths
subsamples_dir = os.path.join(base_path, subsets_filtered_dir_rel)
output_directory = os.path.join(base_path, umap_qt_3d_dir_rel)

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
    scaler = QuantileTransformer(random_state=42)
    features_scaled = scaler.fit_transform(features_imputed)
    
    # Run UMAP with 3 components
    umap_model = umap.UMAP(n_components=3, random_state=42)
    umap_embedding = umap_model.fit_transform(features_scaled)
    
    # Prepare the plot data
    plot_data = pd.DataFrame(umap_embedding, columns=['UMAP-1', 'UMAP-2', 'UMAP-3'])
    plot_data['Concentration'] = metadata['Concentration'].values
    
    # Determine the base file name without extension and directory
    base_file_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # Create a color dictionary
    concentrations = sorted(plot_data['Concentration'].unique())
    color_dict = assign_colors(concentrations)
    
    # Function to create 2D plots (updated to match umap_QT_subsamples_concentrations.py)
    def create_2d_plot(x, y, title, filename):
        plt.figure(figsize=(10, 8))
        sns.scatterplot(
            x=x,
            y=y,
            hue='Concentration',
            data=plot_data,
            palette=color_dict,
            s=9,
            hue_order=sorted(concentrations, reverse=True),
            #legend='brief'
        )
        plt.title(title)
        plt.xlabel(x)
        plt.ylabel(y)
        #plt.legend(title='Concentration', loc='best', bbox_to_anchor=None)
        plt.gca().get_legend().remove()
        plt.savefig(os.path.join(output_directory, filename), dpi=300)
        plt.close()
    
    # Create 2D plots
    create_2d_plot('UMAP-1', 'UMAP-2', 'UMAP plot colored by Concentration', f"umap_{base_file_name}_concentration_dim1_dim2.png")
    create_2d_plot('UMAP-1', 'UMAP-3', 'UMAP plot colored by Concentration', f"umap_{base_file_name}_concentration_dim1_dim3.png")
    create_2d_plot('UMAP-2', 'UMAP-3', 'UMAP plot colored by Concentration', f"umap_{base_file_name}_concentration_dim2_dim3.png")
    
    # Create 3D plot (unchanged)
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    for concentration in concentrations:
        subset = plot_data[plot_data['Concentration'] == concentration]
        ax.scatter(
            subset['UMAP-1'],
            subset['UMAP-2'],
            subset['UMAP-3'],
            c=color_dict[concentration],
            s=7,
            #label=str(concentration)
        )
    
    ax.set_xlabel('UMAP-1')
    ax.set_ylabel('UMAP-2')
    ax.set_zlabel('UMAP-3')
    ax.set_title('3D UMAP plot colored by Concentration')
    #ax.legend(title='Concentration', loc='best', bbox_to_anchor=(1.05, 1), ncol=1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, f"umap_{base_file_name}_concentration_3d.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save the UMAP embedding for later use
    embedding_save_path = os.path.join(output_directory, f"{base_file_name}_umap_embedding_QuantileTransformer_3d.npy")
    np.save(embedding_save_path, umap_embedding)

print("UMAP plots with 3 dimensions (including 3D plot) and embeddings created for all subsamples with QuantileTransformer scaling!")