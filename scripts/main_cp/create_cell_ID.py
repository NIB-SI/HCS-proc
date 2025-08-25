import pandas as pd
import os

def create_unique_id(df, id_name, image_col, object_col):
    df[id_name] = df[image_col].astype(str) + '_' + df[object_col].astype(str)
    return df

# Define the base directory
base_dir = '/PATH/TO/main/results'

files = {
    "Cells": "mod_NEWfn_Cells.txt",
    "partialCells": "mod_NEWfn_partialCells.txt",
    "FilteredNuclei": "mod_NEWfn_FilteredNuclei.txt",
    "RelateMitoCell": "mod_NEWfn_RelateMitoCell.txt",
    "PunctaMito": "mod_NEWfn_PunctaMito.txt",
    "RelateLysoCell": "mod_NEWfn_RelateLysoCell.txt",
    "PunctaLyso": "mod_NEWfn_PunctaLyso.txt",
    "Cytoplasm": "mod_NEWfn_Cytoplasm.txt"
}

files_to_delete = [
    "mod_NEWfn_PunctaLyso.txt",
    "mod_NEWfn_PunctaMito.txt",
    "mod_NEWfn_partialCells.txt"
]

dfs = {key: pd.read_csv(os.path.join(base_dir, files[key]), sep='\t') for key in files}

# Create unique IDs and link columns
# Path 1
dfs['Cells'] = create_unique_id(dfs['Cells'], 'cell_ID', 'ImageNumber', 'ObjectNumber')
dfs['Cells']['link_Parent_partialCells'] = dfs['Cells']['ImageNumber'].astype(str) + '_' + dfs['Cells']['Parent_partialCells'].astype(str)

dfs['partialCells'] = create_unique_id(dfs['partialCells'], 'partialCells_ID', 'ImageNumber', 'ObjectNumber')
dfs['partialCells'] = dfs['partialCells'].merge(dfs['Cells'][['link_Parent_partialCells', 'cell_ID']], left_on='partialCells_ID', right_on='link_Parent_partialCells', how='left')
dfs['partialCells']['link_Parent_FilteredNuclei'] = dfs['partialCells']['ImageNumber'].astype(str) + '_' + dfs['partialCells']['Parent_FilteredNuclei'].astype(str)

dfs['FilteredNuclei'] = create_unique_id(dfs['FilteredNuclei'], 'FilteredNuclei_ID', 'ImageNumber', 'ObjectNumber')
dfs['FilteredNuclei'] = dfs['FilteredNuclei'].merge(dfs['partialCells'][['cell_ID', 'link_Parent_FilteredNuclei']], left_on='FilteredNuclei_ID', right_on='link_Parent_FilteredNuclei', how='left')

# Path 2 
for entity in ['Mito', 'Lyso']:
    relate_entity = f'Relate{entity}Cell'
    puncta_entity = f'Puncta{entity}'

    dfs[relate_entity] = create_unique_id(dfs[relate_entity], f'{relate_entity}_ID', 'ImageNumber', 'ObjectNumber')
    dfs[relate_entity][f'link_Parent_{puncta_entity}'] = dfs[relate_entity]['ImageNumber'].astype(str) + '_' + dfs[relate_entity][f'Parent_{puncta_entity}'].astype(str)

    dfs[puncta_entity] = create_unique_id(dfs[puncta_entity], f'{puncta_entity}_ID', 'ImageNumber', 'ObjectNumber')
    dfs[puncta_entity] = dfs[puncta_entity].merge(dfs[relate_entity][[f'{relate_entity}_ID', f'link_Parent_{puncta_entity}']], left_on=f'{puncta_entity}_ID', right_on=f'link_Parent_{puncta_entity}', how='left')
    dfs[puncta_entity]['cell_ID'] = dfs[puncta_entity]['ImageNumber'].astype(str) + '_' + dfs[puncta_entity]['Parent_Cells'].astype(str)

    dfs[relate_entity] = dfs[relate_entity].merge(dfs[puncta_entity][['cell_ID', f'{relate_entity}_ID']], on=f'{relate_entity}_ID', how='left')

# Path 3
dfs['Cytoplasm']['cell_ID'] = dfs['Cytoplasm']['ImageNumber'].astype(str) + '_' + dfs['Cytoplasm']['Parent_Cells'].astype(str)

# Removing control values at B11 site 5_*
site_pattern = r'5_\d+'

# Convert the 'Metadata_Site' column to string type
for key, df in dfs.items():
    df['Metadata_Site'] = df['Metadata_Site'].astype(str)
    df = df[~((df['Metadata_Well'] == 'B11') & df['Metadata_Site'].str.match(site_pattern))]
    # Save the updated dataframe back to its file
    df.to_csv(os.path.join(base_dir, files[key]), sep='\t', index=False)
    print(f'Updated {files[key]} with new identifiers.')

# Remove the specified files
for filename in files_to_delete:
    file_path = os.path.join(base_dir, filename)
    if os.path.exists(file_path):
        os.remove(file_path)
        print(f'Removed file: {filename}')
    else:
        print(f'File not found and could not be deleted: {filename}')
