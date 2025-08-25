import pandas as pd
import os
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_subsets']['base_path']
final_features_file_rel = config['dim_reduction_subsets']['final_features_file']
subsets_all_days_fs_dir_rel = config['dim_reduction_subsets']['subsets_all_days_fs_dir']

# Construct full paths
file_path = os.path.join(base_path, final_features_file_rel)
output_directory = os.path.join(base_path, subsets_all_days_fs_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Load the data
data = pd.read_csv(file_path, sep='\t')

# Ensure 'Concentration' is of string type
data['Concentration'] = data['Concentration'].astype(str)

# Define the days to be included in subsamples
days_to_include = ['D1', 'D5', 'D7', 'D9']

# Function to create and save subsamples
def create_subsamples(data, days, concentrations, exclude_concentration=None):
    for day in days:
        # Filter out rows for the current day
        day_data = data.loc[data['Metadata_Day'].str.upper() == day].copy()
        print(f"Day {day} data shape after filtering: {day_data.shape}")  # Debug statement
        
        # Ensure 'Concentration' is of string type
        day_data['Concentration'] = day_data['Concentration'].astype(str)
        
        # Initialize a DataFrame to hold the subsample
        subsampled_data = pd.DataFrame()
        
        for concentration in concentrations:
            # Skip the excluded concentration
            if str(concentration) == str(exclude_concentration):
                continue
            
            # Filter data for the current concentration
            concentration_data = day_data.loc[day_data['Concentration'].str.strip() == str(concentration)]
            if concentration_data.shape[0] == 0:
                print(f"Skipping day {day}, concentration {concentration} due to no data")
                continue  # Skip if no data for this concentration
            
            print(f"Day {day}, Concentration {concentration} data shape: {concentration_data.shape}")  # Debug statement
            
            # Sample up to 5000 rows from each concentration group
            subsample = concentration_data.sample(n=min(5000, len(concentration_data)), random_state=1)
            
            # Append to the subsampled_data DataFrame
            subsampled_data = pd.concat([subsampled_data, subsample], ignore_index=True)
        
        # Shuffle the subsampled data
        if not subsampled_data.empty:
            subsampled_data = subsampled_data.sample(frac=1, random_state=1).reset_index(drop=True)
        
        # Define the output file path
        if exclude_concentration:
            new_filename = f"subsample_5k_{day}_noC{exclude_concentration}.txt"
        else:
            new_filename = f"subsample_5k_{day}.txt"
        output_file_path = os.path.join(output_directory, new_filename)
        
        # Save the subsampled data to a new tab-delimited file
        subsampled_data.to_csv(output_file_path, sep='\t', index=False)
        print(f"Subsampled data for day {day} saved to {output_file_path}")  # Debug statement

# Define all concentrations
all_concentrations = list(range(2, 12))  # Assuming concentrations 2 through 11 are inclusive

# Create subsamples including all concentrations
create_subsamples(data, days_to_include, all_concentrations)

# Create subsamples excluding concentration "2"
create_subsamples(data, days_to_include, all_concentrations, exclude_concentration=2)

print("Subsamples created successfully!")