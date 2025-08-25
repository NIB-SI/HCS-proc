import glob
import os

import pandas as pd

# Define the directories using raw strings
source_directory = "/PATH/TO/preprocessing/qc"
destination_directory = "/PATH/TO/preprocessing/qc/remove"

# Read the "imageID.csv" file
imageID_file_path = os.path.join(source_directory, "imageID.csv")
imageID_df = pd.read_csv(
    imageID_file_path,
    usecols=[
        "ImageNumber",
        "Image_FileName_Lys",
        "Image_FileName_Mito",
        "Image_FileName_Nuclei",
        "Image_Metadata_Biorep",
        "Image_Metadata_Chemical",
        "Image_Metadata_Day",
        "Image_Metadata_Site",
        "Image_Metadata_Well",
        "Image_PathName_Lys",
        "Image_PathName_Mito",
        "Image_PathName_Nuclei",
    ],
)

# Find all the "tagged_" files in the directory
tagged_files = glob.glob(os.path.join(source_directory, "tagged_*.csv"))

# Initialize a list to hold dataframes
merged_dataframes = []

# Loop through the tagged files and merge them with the imageID_df
for file in tagged_files:
    tagged_df = pd.read_csv(file)

    # Clean the column names by replacing new lines and trimming spaces
    tagged_df.columns = [col.replace("\n", " ").strip() for col in tagged_df.columns]

    # Columns to drop, only if they exist in the DataFrame
    columns_to_drop = [
        "Image_Metadata_Biorep",
        "Image_Metadata_Well",
        "Total Image Count",
    ]
    columns_to_drop = [col for col in columns_to_drop if col in tagged_df.columns]

    # Merge while excluding unwanted columns
    merged_df = pd.merge(
        imageID_df,
        tagged_df.drop(columns=columns_to_drop),
        on="ImageNumber",
        how="left",
    )
    merged_dataframes.append(merged_df)

# Concatenate all merged dataframes into a single dataframe
flag_image_list_df = pd.concat(merged_dataframes, ignore_index=True)

# Clean the column names by removing quotes and replacing new lines and trimming spaces
flag_image_list_df.columns = [
    col.replace('"', "").replace("\n", " ").strip()
    for col in flag_image_list_df.columns
]

# Save the cleaned "flag_image_list.csv"
flag_image_list_path = os.path.join(destination_directory, "flag_image_list.csv")
flag_image_list_df.to_csv(flag_image_list_path, index=False)

# Adjust the training_columns to match the actual column names in flag_image_list_df
training_columns = [
    "Artefacts  Image Count",
    "Empty  Image Count",
    "Oversaturated  Image Count"
]

# Ensure the column names match exactly, including spaces
training_columns = [
    col for col in training_columns if col in flag_image_list_df.columns
]

# Proceed only if all the required columns are present
if len(training_columns) == 3:
    mask = (flag_image_list_df[training_columns] == 1).any(axis=1)
    training_df = flag_image_list_df.loc[
        mask,
        [
            "Image_FileName_Lys",
            "Image_FileName_Mito",
            "Image_FileName_Nuclei",
            "Image_PathName_Lys",
        ],
    ]

    # Combine 'Image_FileName_' columns into one 'remove' column and remove duplicates
    remove_df = pd.DataFrame(
        training_df[
            ["Image_FileName_Lys", "Image_FileName_Mito", "Image_FileName_Nuclei"]
        ].values.flatten(),
        columns=["remove"],
    )
    remove_df = remove_df.drop_duplicates().reset_index(drop=True)

    # Save this as 'remove_image_list.csv'
    remove_file_path = os.path.join(destination_directory, "remove_image_list.csv")
    remove_df.to_csv(remove_file_path, index=False)

    print(f"Remove Image List saved to {remove_file_path}")
else:
    print("Not all required training columns were found in the DataFrame.")
