import pandas as pd
import numpy as np

file_path = 'PATH/TO/main/results/cell_ID_pooled_median.txt'
data = pd.read_csv(file_path, sep='\t')

# Define metadata columns
metadata_columns = [
    'Concentration', 'counts_Cells', 'counts_Cytoplasm', 'counts_FilteredNuclei',
    'counts_RelateLysoCell', 'counts_RelateMitoCell', 'Metadata_Well', 'Metadata_Day',
    'Metadata_Biorep', 'Tech_replica', 'Day_Well_BR'
]

# Confirm that the expected metadata columns are present
missing_columns = [col for col in metadata_columns if col not in data.columns]
if missing_columns:
    raise ValueError(f"Missing columns in the dataset: {missing_columns}")

# Identify numeric features by excluding metadata columns and selecting only numeric types
features = data.drop(columns=metadata_columns).select_dtypes(include=[np.number]).columns.tolist()

# Function to manually calculate the median and MAD from a Series
def calculate_manual_mad(x):
    median = x.median()
    deviations = np.abs(x - median)
    mad = deviations.median()
    return median, mad

# Process each feature for robust Z-score and row normalization
for feature in features:
    # Group by Metadata_Day and Metadata_Biorep, apply manual MAD calculation
    results = data.groupby(['Metadata_Day', 'Metadata_Biorep'])[feature].apply(calculate_manual_mad).apply(pd.Series)
    results.columns = ['Median_' + feature, 'MAD_' + feature]
    data = data.merge(results, on=['Metadata_Day', 'Metadata_Biorep'], how='left')
    data['rZs_' + feature] = 0.6745 * (data[feature] - data['Median_' + feature]) / data['MAD_' + feature]

# Prepare DataFrames for saving: one with robust Z-scores and one with row normalizations
robust_z_scores = data[metadata_columns + [col for col in data.columns if col.startswith('rZs_')]]
#row_normalizations = data[metadata_columns + [col for col in data.columns if col.startswith('RN_')]]

# Save the results to separate files
output_path_rzs = file_path.replace('.txt', '_z_standardization.txt')
robust_z_scores.to_csv(output_path_rzs, sep='\t', index=False)
print('Robust Z-scores saved to:', output_path_rzs)
