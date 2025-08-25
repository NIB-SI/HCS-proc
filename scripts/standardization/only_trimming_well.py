import pandas as pd
import numpy as np
import os
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['standardization']['base_path']
output_dir_rel = config['standardization']['output_dir']
feature_selection_dir_rel = config['standardization']['feature_selection_dir']

# Construct full paths
output_dir = os.path.join(base_path, output_dir_rel)
feature_selection_dir = os.path.join(base_path, feature_selection_dir_rel)

# Define input and output file paths
input_file = os.path.join(output_dir, "cell_ID_pooled_median_row_plate_standardization.txt")
output_file_well = os.path.join(feature_selection_dir, "cell_ID_pooled_median_row_plate_standardization_trimmed_well_2.5_97.5.txt")

# Create output directory if it doesn't exist
os.makedirs(feature_selection_dir, exist_ok=True)

# Set trimming percentage
lower_percentile = 2.5
upper_percentile = 97.5

# Read the data
data = pd.read_csv(input_file, sep="\t")

# Strip spaces and normalize case in the key columns
data["Metadata_Day"] = data["Metadata_Day"].astype(str).str.strip().str.upper()
data["Metadata_Biorep"] = data["Metadata_Biorep"].str.strip().str.upper()
data["Concentration"] = data["Concentration"].astype(str).str.strip()
data["Metadata_Well"] = data["Metadata_Well"].astype(str).str.strip().str.upper()
data["Day_Well_BR"] = data["Day_Well_BR"].astype(str).str.strip().str.upper()
data["cell_ID"] = data["cell_ID"].astype(str).str.strip().str.upper()

# Define the metadata columns
metadata_columns = [
    "Concentration", "counts_Cells", "counts_Cytoplasm", "counts_FilteredNuclei", 
    "Metadata_Well", "Metadata_Day", "Metadata_Biorep", "Tech_replica", "Day_Well_BR", "cell_ID"
]

# Prepare a list to store the results
trimmed_data_list = []

# Function to trim the extreme values
def trim_extremes(values, lower_percentile=0, upper_percentile=100):
    if len(values) == 0:
        return values
    lower_bound = np.nanpercentile(values, lower_percentile)
    upper_bound = np.nanpercentile(values, upper_percentile)
    return values[(values >= lower_bound) & (values <= upper_bound)]

# Trim the data at the well level
grouped = data.groupby(["Metadata_Day", "Metadata_Biorep", "Metadata_Well"])
for (day, biorep, well), group in grouped:
    trimmed_group = group.copy()
    for feature in group.columns.difference(metadata_columns):
        trimmed_values = trim_extremes(group[feature].dropna(), lower_percentile, upper_percentile)
        trimmed_group[feature] = trimmed_values
    trimmed_data_list.append(trimmed_group)

# Concatenate all trimmed data into a single DataFrame
trimmed_data_well = pd.concat(trimmed_data_list)

# Save the trimmed data to a file
trimmed_data_well.to_csv(output_file_well, sep="\t", index=False)

print("Trimming by Well level completed and saved to", output_file_well)