import os
import logging
import configparser
from typing import List, Tuple

import dask.dataframe as dd
import pandas as pd
import numpy as np

# Set up config file
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, 'config_day_well_br.ini')

# Read configuration
config = configparser.ConfigParser()
config.read(config_path)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

BASE_PATH = config['Paths']['base_path']
FILENAMES = config['Files']['filenames'].split(',')
ID_COLUMN = config['Columns']['id_columns']

def process_file(filepath: str, file_identifier: str) -> pd.DataFrame:
    ddf = dd.read_csv(filepath, sep='\t')
    
    if ID_COLUMN not in ddf.columns:
        raise ValueError(f"{ID_COLUMN} column not found in {filepath}")

    numeric_cols = ddf.select_dtypes(include=[np.number]).columns
    ddf = ddf[[ID_COLUMN] + list(numeric_cols)]
    
    # Group by ID_COLUMN and calculate median
    grouped = ddf.groupby(ID_COLUMN)[numeric_cols].median().compute()
    
    # Rename columns to include file identifier
    grouped = grouped.rename(columns={col: f"{col}_{file_identifier}" for col in grouped.columns})
    
    return grouped

def generate_counts(filepath: str, file_identifier: str) -> pd.DataFrame:
    ddf = dd.read_csv(filepath, sep='\t')
    counts = ddf.groupby(ID_COLUMN).size().compute().reset_index(name=f'counts_{file_identifier}')
    return counts

def load_and_process_files(base_path: str, filenames: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    processed_dfs = []
    count_files = []

    for filename in filenames:
        file_identifier = filename.split('_')[2].split('.')[0]
        filepath = os.path.join(base_path, filename)
        
        if not os.path.exists(filepath):
            logging.warning(f"File not found: {filepath}")
            continue

        processed_df = process_file(filepath, file_identifier)
        processed_dfs.append(processed_df)

        # Generate counts
        counts = generate_counts(filepath, file_identifier)
        count_file = os.path.join(base_path, f'intermediate_counts_{file_identifier}.txt')
        counts.to_csv(count_file, sep='\t', index=False)
        count_files.append(count_file)

    if not processed_dfs:
        logging.error("No valid files found to process.")
        return pd.DataFrame(), []

    # Combine all processed dataframes
    combined_df = pd.concat(processed_dfs, axis=1)
    
    return combined_df, count_files

    if not processed_dfs:
        logging.error("No valid files found to process.")
        return pd.DataFrame(), []

    # Combine all processed dataframes
    combined_df = pd.concat(processed_dfs, axis=1)
    
    return combined_df, count_files

def add_counts(data: pd.DataFrame, count_files: List[str]) -> pd.DataFrame:
    for count_file in count_files:
        # Change this line if necessary
        file_identifier = os.path.basename(count_file).split('_')[2].split('.')[0]
        counts = pd.read_csv(count_file, sep='\t')
        data = pd.merge(data, counts, on=ID_COLUMN, how='left')
    return data

def add_metadata(data: pd.DataFrame, base_path: str) -> pd.DataFrame:
    cells_file = [f for f in FILENAMES if 'Cells' in f][0]
    filepath = os.path.join(base_path, cells_file)
    
    metadata = pd.read_csv(filepath, sep='\t', usecols=[ID_COLUMN, 'Metadata_Well', 'Metadata_Day', 'Metadata_Biorep'])
    metadata = metadata.drop_duplicates(ID_COLUMN)
    
    data = pd.merge(data, metadata, on=ID_COLUMN, how='left')
    
    data['Tech_replica'] = data['Metadata_Well'].str.extract(r'([A-Z])')
    data['Concentration'] = pd.to_numeric(data['Metadata_Well'].str.extract(r'(\d+)')[0], errors='coerce')
    
    return data

def clean_columns(data: pd.DataFrame) -> pd.DataFrame:
    columns_to_exclude = ['ImageNumber', 'ObjectNumber', 'Metadata_Channel',
                          'Metadata_Chemical', 'Metadata_Frame', 'Metadata_Series', 
                          'Metadata_Site', 'Subfolder']
    
    columns_to_keep = [ID_COLUMN] + [col for col in data.columns if not (
        col.endswith('_ID') or col.startswith('link_') or col.startswith('Parent_') or 
        col.startswith('FileName_') or col.startswith('PathName_') or 
        col.startswith('Metadata_Channel') or col.startswith('Metadata_Frame') or
        col.startswith('Metadata_Series') or col.startswith('Metadata_Site') or
        col in columns_to_exclude)]
    
    return data[columns_to_keep]

def final_cleanup(data: pd.DataFrame) -> pd.DataFrame:
    # Remove ObjectNumber_* and ImageNumber_* columns
    columns_to_keep = [col for col in data.columns if not (col.startswith('ObjectNumber_') or col.startswith('ImageNumber_'))]
    data = data[columns_to_keep]
    
    # Handle Day_Well_BR columns
    day_well_br_columns = [col for col in data.columns if col == 'Day_Well_BR']
    if len(day_well_br_columns) > 1:
        # Find the index of the first 'Day_Well_BR' column
        first_day_well_br_index = data.columns.get_loc('Day_Well_BR')
        if isinstance(first_day_well_br_index, (int, np.integer)):
            # Keep only the first Day_Well_BR column
            columns_to_keep = [col for col in data.columns if col != 'Day_Well_BR'] + ['Day_Well_BR']
            data = data[columns_to_keep]
            # Ensure 'Day_Well_BR' is at its original position
            cols = data.columns.tolist()
            cols.remove('Day_Well_BR')
            cols.insert(first_day_well_br_index, 'Day_Well_BR')
            data = data[cols]
        else:
            # If there are multiple matches, just keep the first one
            data = data.loc[:, ~data.columns.duplicated(keep='first')]
    
    return data

def main():
    try:
        logging.info(f"Processing data for {ID_COLUMN}")
        
        pooled_data, count_files = load_and_process_files(BASE_PATH, FILENAMES)
        logging.info(f"Pooled data shape: {pooled_data.shape}")
        
        pooled_data_with_counts = add_counts(pooled_data, count_files)
        logging.info(f"Data shape after adding counts: {pooled_data_with_counts.shape}")
        
        data_with_metadata = add_metadata(pooled_data_with_counts, BASE_PATH)
        logging.info(f"Data shape after adding metadata: {data_with_metadata.shape}")
        
        cleaned_data = clean_columns(data_with_metadata)
        logging.info(f"Data shape after cleaning columns: {cleaned_data.shape}")
        
        final_data = final_cleanup(cleaned_data)
        logging.info(f"Final data shape after removing ObjectNumber, ImageNumber, and duplicate Day_Well_BR: {final_data.shape}")

        output_file = os.path.join(BASE_PATH, f'{ID_COLUMN}_pooled_median.txt')
        final_data.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Output saved to {output_file}")

        # Clean up intermediate count files
        for file in count_files:
            os.remove(file)
        logging.info("Intermediate count files cleaned up")

    except Exception as e:
        logging.error(f"Error in main function: {str(e)}")
        logging.error(f"Error occurred on line: {e.__traceback__.tb_lineno}")
        raise

if __name__ == "__main__":
    main()