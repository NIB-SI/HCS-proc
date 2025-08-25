import pandas as pd
import os
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_subsets']['base_path']
final_features_file_rel = config['dim_reduction_subsets']['final_features_file']
subsets_days_per_conc_dir_rel = config['dim_reduction_subsets']['subsets_days_per_conc_dir']

# Construct full paths
file_path = os.path.join(base_path, final_features_file_rel)
output_directory = os.path.join(base_path, subsets_days_per_conc_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Load the data
data = pd.read_csv(file_path, sep='\t')

# Ensure 'Concentration' is of string type
data['Concentration'] = data['Concentration'].astype(str)

# Define the days to be included in subsamples
days_to_include = ['D1', 'D5', 'D7', 'D9']

# Define concentration-specific cutoffs
def get_cutoff_for_concentration(concentration):
    if concentration == 2:
        return 5000
    elif concentration == 3:
        return 7000
    elif concentration in range(4, 12):  # 4-11
        return 10000
    else:
        return 10000  # default for any other concentration

# Function to create and save subsamples
def create_subsamples(data, days, concentrations):
    for concentration in concentrations:
        # Get the appropriate cutoff for this concentration
        cutoff = get_cutoff_for_concentration(concentration)
        
        # Filter out rows for the current concentration
        concentration_data = data.loc[data['Concentration'].str.strip() == str(concentration)].copy()
        print(f"Concentration {concentration} data shape after filtering: {concentration_data.shape}")  # Debug statement
        
        # Initialize a DataFrame to hold the subsample
        subsampled_data = pd.DataFrame()
        
        for day in days:
            # Filter data for the current day
            day_data = concentration_data.loc[concentration_data['Metadata_Day'].str.upper() == day]
            if day_data.shape[0] == 0:
                print(f"Skipping concentration {concentration}, day {day} due to no data")
                continue  # Skip if no data for this day
            
            print(f"Concentration {concentration}, Day {day} data shape: {day_data.shape}")  # Debug statement
            
            # Sample up to the specified cutoff from each day group
            subsample = day_data.sample(n=min(cutoff, len(day_data)), random_state=1)
            
            # Append to the subsampled_data DataFrame
            subsampled_data = pd.concat([subsampled_data, subsample], ignore_index=True)
        
        # Shuffle the subsampled data
        if not subsampled_data.empty:
            subsampled_data = subsampled_data.sample(frac=1, random_state=1).reset_index(drop=True)
        
        # Define the output file path
        new_filename = f"subsample_5k_C{concentration}.txt"
        output_file_path = os.path.join(output_directory, new_filename)
        
        # Save the subsampled data to a new tab-delimited file
        subsampled_data.to_csv(output_file_path, sep='\t', index=False)
        print(f"Subsampled data for concentration {concentration} saved to {output_file_path}")  # Debug statement

# Define all concentrations
all_concentrations = list(range(2, 12))  # Assuming concentrations 2 through 11 are inclusive

# Create subsamples for each concentration
create_subsamples(data, days_to_include, all_concentrations)

print("Subsamples created successfully!")