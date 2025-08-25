import pandas as pd
import os
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_subsets']['base_path']
final_features_file_rel = config['dim_reduction_subsets']['final_features_file']
counts_output_dir_rel = config['dim_reduction_subsets']['counts_output_dir']

# Construct full paths
file_path = os.path.join(base_path, final_features_file_rel)
output_dir = os.path.join(base_path, counts_output_dir_rel)

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load the data
data = pd.read_csv(file_path, sep='\t')

# 1. Count for each unique Concentration per Metadata_Day
day_concentration_counts = data.groupby(['Metadata_Day', 'Concentration']).size().reset_index(name='Counts')
print("Counts for each unique Concentration per Metadata_Day:")
print(day_concentration_counts)
day_concentration_counts.to_csv(os.path.join(output_dir, 'day_concentration_counts.txt'), sep='\t', index=False)

# 2. Count for each Metadata_Biorep per Metadata_Day
day_biorep_counts = data.groupby(['Metadata_Day', 'Metadata_Biorep']).size().reset_index(name='Counts')
print("\nCounts for each Metadata_Biorep per Metadata_Day:")
print(day_biorep_counts)
day_biorep_counts.to_csv(os.path.join(output_dir, 'day_biorep_counts.txt'), sep='\t', index=False)

# 3. Count for each Metadata_Biorep per Concentration
concentration_biorep_counts = data.groupby(['Concentration', 'Metadata_Biorep']).size().reset_index(name='Counts')
print("\nCounts for each Metadata_Biorep per Concentration:")
print(concentration_biorep_counts)
concentration_biorep_counts.to_csv(os.path.join(output_dir, 'concentration_biorep_counts.txt'), sep='\t', index=False)

# 4. Count for each Metadata_Biorep per Metadata_Day and Concentration
day_concentration_biorep_counts = data.groupby(['Metadata_Day', 'Concentration', 'Metadata_Biorep']).size().reset_index(name='Counts')
print("\nCounts for each Metadata_Biorep per Metadata_Day and Concentration:")
print(day_concentration_biorep_counts)
day_concentration_biorep_counts.to_csv(os.path.join(output_dir, 'day_concentration_biorep_counts.txt'), sep='\t', index=False)

# 5. Count for each Metadata_Day
day_counts = data['Metadata_Day'].value_counts().reset_index()
day_counts.columns = ['Metadata_Day', 'Counts']
day_counts = day_counts.sort_values('Metadata_Day')
print("\nCounts for each Metadata_Day:")
print(day_counts)
day_counts.to_csv(os.path.join(output_dir, 'day_counts.txt'), sep='\t', index=False)

# 6. Count for each Concentration
concentration_counts = data['Concentration'].value_counts().reset_index()
concentration_counts.columns = ['Concentration', 'Counts']
concentration_counts = concentration_counts.sort_values('Concentration')
print("\nCounts for each Concentration:")
print(concentration_counts)
concentration_counts.to_csv(os.path.join(output_dir, 'concentration_counts.txt'), sep='\t', index=False)

print("\nAll count files have been saved in the directory:", output_dir)