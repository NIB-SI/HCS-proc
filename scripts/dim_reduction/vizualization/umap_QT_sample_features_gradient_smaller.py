import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from sklearn.impute import SimpleImputer
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from PIL import Image, ImageOps
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_visualization']['base_path']
sample_features_file_rel = config['dim_reduction_visualization']['sample_features_file']
sample_embedding_file_rel = config['dim_reduction_visualization']['sample_embedding_file']

# Construct full paths
file_path = os.path.join(base_path, sample_features_file_rel)
embedding_file_path = os.path.join(base_path, sample_embedding_file_rel)

# Function to normalize a series between its 5th and 95th percentile
def normalize_series(s):
    # Convert to numeric, coercing errors to NaN
    s = pd.to_numeric(s, errors='coerce')
    
    # Drop NaN values that were non-numeric before normalization
    s = s.dropna()

    lower_bound = s.quantile(0.05)
    upper_bound = s.quantile(0.95)
    s = s.clip(lower_bound, upper_bound)
    s = (s - lower_bound) / (upper_bound - lower_bound)
    return s, lower_bound, upper_bound

# Load the data
data = pd.read_csv(file_path, sep='\t')

# Get base file name for output
base_file_name = os.path.splitext(os.path.basename(file_path))[0]

# Create output directory for plots
plot_directory = os.path.join(base_path, f"dim_reduction/results/umap/QT/plots_{base_file_name}_features")
os.makedirs(plot_directory, exist_ok=True)

# Randomize the order of rows in the dataset
data = data.sample(frac=1, random_state=42).reset_index(drop=True)

# Load the saved UMAP embedding (which was generated after the data was shuffled with random_state=42)
umap_embedding = np.load(embedding_file_path)

# Separate features and metadata
metadata_columns = ['Concentration', 'counts_Cells', 'counts_Cytoplasm', 'counts_FilteredNuclei',
                    'Metadata_Well', 'Metadata_Day', 'Metadata_Biorep', 'Tech_replica', 'Day_Well_BR', 'cell_ID']

# Define features to exclude
#features_to_exclude = [
#    'rp_norm_Mean_PunctaLyso_Distance_Minimum_Cytoplasm_FilteredNuclei',
#    'rp_norm_AreaShape_FormFactor_RelateMitoCell',
#    'rp_norm_Number_Object_Number_Cells',
#    'rp_norm_AreaShape_Compactness_RelateLysoCell',
#    'rp_norm_Texture_AngularSecondMoment_GrayLys_3_00_256_RelateLysoCell',
#    'rp_norm_Mean_PunctaLyso_Number_Object_Number_Cells'
#]

# Drop metadata columns and excluded features
#features = data.drop(columns=metadata_columns + features_to_exclude)
features = data.drop(columns=metadata_columns)

# Replace infinite values with NaN
features = features.replace([np.inf, -np.inf], np.nan)

# Function to clean up UMAP axes
def setup_clean_umap_axes(ax):
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # Remove the spines (borders)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return ax

# Function to clean up density plot axes but keep x-axis
def setup_clean_density_axes(ax):
    ax.set_ylabel('')
    ax.set_yticks([])
    ax.set_yticklabels([])
    # Keep x-axis ticks and labels
    
    # Only keep bottom spine
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_linewidth(0.5)  # Make it thin
    
    return ax

# Create a plot for each feature
for column in features.columns:
    if pd.api.types.is_numeric_dtype(data[column]):
        print(f"Processing feature: {column}")
        
        # Step 1: Prepare the UMAP coordinates DataFrame
        plot_data = pd.DataFrame({
            'UMAP-1': umap_embedding[:, 0],
            'UMAP-2': umap_embedding[:, 1]
        })
        
        # Step 2: Normalize feature data and handle missing values
        normalized_data, lower_bound, upper_bound = normalize_series(data[column])
        
        # Create a mask for missing values
        plot_data[column] = np.nan
        valid_indices = normalized_data.index
        common_indices = [i for i in valid_indices if i < len(plot_data)]
        for i in common_indices:
            plot_data.loc[i, column] = normalized_data.loc[i]
        
        # Get valid data for distribution plot
        valid_data = data[column].replace([np.inf, -np.inf], np.nan).dropna()
        valid_data_clipped = valid_data.clip(lower_bound, upper_bound)
        
        missing_mask = plot_data[column].isna()
        non_missing_data = plot_data[~missing_mask] if missing_mask.any() else plot_data
        
        # ----- Generate UMAP plot -----
        umap_fig = plt.figure(figsize=(10, 8))
        ax_umap = umap_fig.add_subplot(111)
        
        # Get the exact plot limits from the first figure to maintain consistent sizing
        if missing_mask.any():
            sns.scatterplot(
                x='UMAP-1',
                y='UMAP-2',
                data=plot_data[missing_mask],
                color='#E0E0E0',
                s=9,
                alpha=1.0,
                edgecolor=None,
                ax=ax_umap
            )
        
        if len(non_missing_data) > 0:
            sns.scatterplot(
                x='UMAP-1',
                y='UMAP-2',
                hue=column,
                data=non_missing_data,
                palette='inferno',
                s=9,
                alpha=1.0,
                legend=False,
                ax=ax_umap
            )
        
        # Store limits before removing axes
        x_min, x_max = ax_umap.get_xlim()
        y_min, y_max = ax_umap.get_ylim()
        
        # Clean up the axes - remove labels, ticks, etc.
        setup_clean_umap_axes(ax_umap)
        
        # Title at the top with minimal padding
        plt.title(f'UMAP plot colored by {column}', pad=10)
        
        # Make sure to keep the exact same limits after cleanup
        ax_umap.set_xlim(x_min, x_max)
        ax_umap.set_ylim(y_min, y_max)
        
        # Save temporary UMAP plot with minimum padding
        temp_umap_path = os.path.join(plot_directory, "temp_umap.png")
        plt.savefig(temp_umap_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()
        
        # ----- Generate density plot -----
        density_fig = plt.figure(figsize=(10, 1.5))
        ax_dist = density_fig.add_subplot(111)
        
        # Format x-ticks with fewer decimal places if the range is small
        value_range = upper_bound - lower_bound
        if abs(value_range) < 0.1:
            decimal_places = 3
        elif abs(value_range) < 1:
            decimal_places = 2
        else:
            decimal_places = 1
        
        # Get the KDE data for manual coloring
        try:
            kde = sns.kdeplot(
                x=valid_data_clipped,
                ax=ax_dist,
                color='black',
                fill=False,
                linewidth=1
            )
            
            # Get the line data from the KDE plot
            if len(kde.get_lines()) > 0:
                line = kde.get_lines()[-1]
                x_kde = line.get_xdata()
                y_kde = line.get_ydata()
            else:
                print(f"Warning: No KDE lines generated for {column}. This feature may have zero variance. Skipping density plot.")
                # Create dummy data for a flat line at y=0
                x_kde = np.array([lower_bound, upper_bound])
                y_kde = np.array([0, 0])
        except Exception as e:
            print(f"Warning: Error generating KDE for {column}: {str(e)}. Skipping density plot.")
            # Create dummy data for a flat line at y=0
            x_kde = np.array([lower_bound, upper_bound])
            y_kde = np.array([0, 0])
        
        # Create a colormap mapping from data values to normalized positions (0-1)
        norm = plt.Normalize(lower_bound, upper_bound)
        cmap = plt.cm.inferno
        x_norm = (x_kde - lower_bound) / (upper_bound - lower_bound)
        
        # Create overlapping segments with antialiasing to prevent visible borders
        for i in range(len(x_kde) - 1):
            # For each segment, create a polygon with overlap
            overlap = 0.1  # 10% overlap between segments
            
            if i < len(x_kde) - 2:  # Not the last segment
                next_x = x_kde[i+1]
                segment_width = next_x - x_kde[i]
                end_x = next_x + (overlap * segment_width)
            else:  # Last segment
                end_x = x_kde[i+1]
            
            if i > 0:  # Not the first segment
                prev_x = x_kde[i-1]
                segment_width = x_kde[i] - prev_x
                start_x = x_kde[i] - (overlap * segment_width)
            else:  # First segment
                start_x = x_kde[i]
            
            segment_x = [start_x, end_x, end_x, start_x]
            segment_y = [0, 0, y_kde[i+1], y_kde[i]]
            
            # Get color based on x position
            color1 = cmap(x_norm[i])
            color2 = cmap(x_norm[i+1])
            avg_color = 0.5 * np.array(color1) + 0.5 * np.array(color2)
            
            # Plot the polygon with antialiased edges and no edge color
            ax_dist.fill(segment_x, segment_y, color=avg_color, antialiased=True, linewidth=0, edgecolor='none')
        
        # Redraw the black line on top for clarity
        ax_dist.plot(x_kde, y_kde, 'k-', linewidth=1)
        
        # Format x-axis ticks to use consistent decimal places
        ax_dist.xaxis.set_major_formatter(plt.FormatStrFormatter(f'%.{decimal_places}f'))
        
        # Clean up the density plot axes but KEEP the x-axis
        setup_clean_density_axes(ax_dist)
        
        # Add feature value label at the bottom with minimal padding
        ax_dist.set_xlabel(f'Feature value ({column})', labelpad=5)
        
        # Save temporary density plot with minimal padding
        temp_density_path = os.path.join(plot_directory, "temp_density.png")
        plt.savefig(temp_density_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()
        
        # ----- Combine the images using PIL -----
        umap_img = Image.open(temp_umap_path)
        density_img = Image.open(temp_density_path)
        
        # Resize density image to exactly match width of UMAP image
        density_img = density_img.resize((umap_img.width, density_img.height), Image.LANCZOS)
        
        # Create a new blank image to hold both
        combined_height = umap_img.height + density_img.height
        combined_img = Image.new('RGB', (umap_img.width, combined_height), (255, 255, 255))
        
        # Paste the UMAP image at the top
        combined_img.paste(umap_img, (0, 0))
        
        # Paste the density image at the bottom
        combined_img.paste(density_img, (0, umap_img.height))
        
        # Add a thin black line between the two plots
        for x in range(umap_img.width):
            combined_img.putpixel((x, umap_img.height), (0, 0, 0))
        
        # Save the combined image
        combined_path = os.path.join(plot_directory, f"umap_{base_file_name}_{column}.png")
        combined_img.save(combined_path, "PNG")
        
        # Clean up temporary files
        os.remove(temp_umap_path)
        os.remove(temp_density_path)

print(f"All plots have been saved to {plot_directory}")