import os
import pandas as pd

# Define the path to the main results directory
main_dir = "/PATH/TO/main/results"

# File names of interest for combining
file_names = ["NEWfn_Cells.txt", "NEWfn_Cytoplasm.txt", "NEWfn_FilteredNuclei.txt",
              "NEWfn_RelateLysoCell.txt", "NEWfn_RelateMitoCell.txt", "NEWfn_partialCells.txt", 
			  "NEWfn_PunctaLyso.txt", "NEWfn_PunctaMito.txt"]

# File-specific strings to look for in column names for filtering
file_strings = {
    'combined_NEWfn_FilteredNuclei.txt': ["Children", "Location", "GrayLys", "GrayMito", "_01_", "_02_", "_03_"],
    'combined_NEWfn_RelateLysoCell.txt': ["Children", "Location", "GrayNuclei", "GrayMito", "_01_", "_02_", "_03_"],
    'combined_NEWfn_PunctaLyso.txt.txt': ["Children", "Location", "GrayNuclei", "GrayMito", "_01_", "_02_", "_03_"],
    'combined_NEWfn_PunctaMito.txt': ["Children", "Location", "GrayNuclei", "GrayLys", "_01_", "_02_", "_03_"],
	'combined_NEWfn_RelateMitoCell.txt': ["Children", "Location", "GrayNuclei", "GrayLys", "_01_", "_02_", "_03_"],
    'combined_NEWfn_Cells.txt': ["Children", "Location", "_01_", "_02_", "_03_"],
    'combined_NEWfn_Cytoplasm.txt': ["Children", "Location", "_01_", "_02_", "_03_"],
}

def filter_and_save_combined_data(combined_data, file_name):
    original_file_name = 'combined_' + file_name
    strings_to_remove = file_strings.get(original_file_name, [])
    columns_to_drop = [col for col in combined_data.columns if any(substring in col for substring in strings_to_remove)]
    combined_data.drop(columns=columns_to_drop, inplace=True)
    
    # Save the modified DataFrame with 'mod_' prefix
    combined_data.to_csv(os.path.join(main_dir, f'mod_{file_name}'), sep='\t', index=False)

for file_name in file_names:
    combined_data = pd.DataFrame()
    for subfolder in os.listdir(main_dir):
        subfolder_path = os.path.join(main_dir, subfolder)
        if os.path.isdir(subfolder_path) and file_name in os.listdir(subfolder_path):
            file_path = os.path.join(subfolder_path, file_name)
            data = pd.read_csv(file_path, sep="\t", header=0)

            # Add the subfolder information as a new column
            data["Subfolder"] = subfolder

            # Append to the combined data
            combined_data = pd.concat([combined_data, data], ignore_index=True)

    # Filter and save the combined data
    filter_and_save_combined_data(combined_data, file_name)

print("All specified files have been processed and saved with a 'mod_' prefix without saving intermediate combined files.")
