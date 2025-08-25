import pandas as pd
import numpy as np
import configparser
import os

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['standardization']['base_path']
input_file_rel = config['standardization']['input_file']
output_dir_rel = config['standardization']['output_dir']

# Construct full paths
input_file = os.path.join(base_path, input_file_rel)
output_dir = os.path.join(base_path, output_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load data from the file
data = pd.read_csv(input_file, sep='\t')

# Define metadata columns
metadata_columns = [
    'Concentration', 'counts_Cells', 'counts_Cytoplasm', 'counts_FilteredNuclei',
    'counts_RelateLysoCell', 'counts_RelateMitoCell', 'Metadata_Well', 'Metadata_Day',
    'Metadata_Biorep', 'Tech_replica', 'Day_Well_BR', 'cell_ID'
]

# Confirm that the expected metadata columns are present
missing_columns = [col for col in metadata_columns if col not in data.columns]
if missing_columns:
    raise ValueError(f"Missing columns in the dataset: {missing_columns}")

# Identify numeric features by excluding metadata columns and selecting only numeric types
features = data.drop(columns=metadata_columns).select_dtypes(include=[np.number]).columns.tolist()

# Process each feature
for feature in features:
    # Calculate Median_row for each group
    data['Median_row_' + feature] = data.groupby(['Metadata_Day', 'Metadata_Biorep', 'Tech_replica'])[feature].transform('median')
    
    # Calculate row normalization
    data['row_normalization_' + feature] = data[feature] - data['Median_row_' + feature]
    
    # Calculate Median_day for each Metadata_Day
    data['Median_day_' + feature] = data.groupby('Metadata_Day')['row_normalization_' + feature].transform('median')
    
    # Calculate row plate normalization
    data['rp_norm_' + feature] = data['row_normalization_' + feature] - data['Median_day_' + feature]

# Prepare DataFrame for saving
row_plate_normalizations = data[metadata_columns + [col for col in data.columns if col.startswith('rp_norm_')]]

# Save the results to a new file (construct output path from config)
output_filename = os.path.basename(input_file).replace('.txt', '_row_plate_standardization.txt')
output_path_rpn = os.path.join(output_dir, output_filename)
row_plate_normalizations.to_csv(output_path_rpn, sep='\t', index=False)
print('Row plate normalizations saved to:', output_path_rpn)