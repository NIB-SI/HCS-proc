import os
import logging
import configparser
from typing import List, Tuple
import shutil

import dask.dataframe as dd
import pandas as pd
import numpy as np

# Set up config file
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, 'config.ini')

# Read configuration
config = configparser.ConfigParser()
config.read(config_path)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

BASE_PATH = config['Paths']['base_path']
FILENAMES = config['Files']['filenames'].split(',')

def process_chunk(chunk: dd.DataFrame, file_identifier: str) -> dd.DataFrame:
    try:
        # Filter only the relevant columns
        columns_to_exclude = ['ImageNumber', 'ObjectNumber', 'Metadata_Biorep', 'Metadata_Channel',
                              'Metadata_Chemical', 'Metadata_Day', 'Metadata_Frame', 'Metadata_Series', 
                              'Metadata_Site', 'Metadata_Well', 'Subfolder', 'Day_Well_BR']
        columns_to_keep = ['cell_ID'] + [col for col in chunk.columns if not (
            col.endswith('_ID') or col.startswith('link_') or col.startswith('Parent_') or 
            col.startswith('FileName_') or col.startswith('PathName_') or col in columns_to_exclude)]

        chunk = chunk[columns_to_keep]

        # Rename the feature columns to include the file identifier, but don't rename the cell_ID
        chunk = chunk.rename(columns={col: f"{col}_{file_identifier}" for col in chunk.columns if col != 'cell_ID'})

        return chunk
    except Exception as e:
        logging.error(f"Error in process_chunk: {str(e)}")
        raise

def load_and_process_files(base_path: str, filenames: List[str], blocksize: str = '100MB') -> Tuple[pd.DataFrame, List[str]]:
    intermediate_files = []
    count_dfs = []

    try:
        for filename in filenames:
            filepath = os.path.join(base_path, filename)
            file_identifier = filename.split('_')[2].split('.')[0]

            if not os.path.exists(filepath):
                logging.warning(f"File not found: {filepath}")
                continue

            ddf = dd.read_csv(filepath, sep='\t', blocksize=blocksize)
            
            processed_chunks = ddf.map_partitions(process_chunk, file_identifier=file_identifier)

            # Calculate counts
            counts = ddf.groupby('cell_ID').size().to_frame(f'counts_{file_identifier}').reset_index()
            count_dfs.append(counts)

            # Write processed chunks to intermediate file
            intermediate_file = os.path.join(base_path, f'cell_ID_intermediate_{file_identifier}.parquet')
            processed_chunks.to_parquet(intermediate_file, write_index=False)
            intermediate_files.append(intermediate_file)

        if not intermediate_files:
            logging.error(f"No valid data found in any of the input files.")
            return pd.DataFrame(), []

        # Read and combine intermediate files
        combined_df = dd.concat([dd.read_parquet(f) for f in intermediate_files])

        # Ensure no duplicate columns
        combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]

        feature_columns = combined_df.select_dtypes(include=[np.number]).columns.difference(['cell_ID'])

        pooled_data = combined_df.groupby('cell_ID')[feature_columns].median().reset_index().compute()

        # Add counts to the pooled data
        for count_df in count_dfs:
            pooled_data = pd.merge(pooled_data, count_df.compute(), on='cell_ID', how='left')

        return pooled_data, intermediate_files
    except Exception as e:
        logging.error(f"Error in load_and_process_files: {str(e)}")
        raise

def add_metadata(base_path: str, pooled_data: pd.DataFrame) -> pd.DataFrame:
    try:
        cells_filepath = os.path.join(base_path, 'mod_NEWfn_Cells.txt')
        cells_df = pd.read_csv(cells_filepath, sep='\t')

        logging.info(f"Cells dataframe shape: {cells_df.shape}")
        logging.info(f"Cells dataframe columns: {cells_df.columns}")
        logging.info(f"Pooled data shape: {pooled_data.shape}")
        logging.info(f"Pooled data columns: {pooled_data.columns}")

        if 'cell_ID' not in cells_df.columns and 'ImageNumber' in cells_df.columns and 'ObjectNumber' in cells_df.columns:
            cells_df['cell_ID'] = cells_df['ImageNumber'].astype(str) + '_' + cells_df['ObjectNumber'].astype(str)

        merge_columns = ['cell_ID', 'Metadata_Well', 'Metadata_Day', 'Metadata_Biorep', 'Day_Well_BR']
        merged_data = pd.merge(pooled_data, cells_df[merge_columns], on='cell_ID', how='left')
        
        logging.info(f"Merged data shape: {merged_data.shape}")
        logging.info(f"Merged data columns: {merged_data.columns}")

        # Extract Tech_replica and Concentration
        if 'Metadata_Well' in merged_data.columns:
            merged_data['Tech_replica'] = merged_data['Metadata_Well'].str.extract(r'([A-Z])')
            merged_data['Concentration'] = pd.to_numeric(merged_data['Metadata_Well'].str.extract(r'(\d+)')[0], errors='coerce')
        else:
            logging.warning("'Metadata_Well' column not found in merged data. Skipping Tech_replica and Concentration extraction.")

        return merged_data
    except Exception as e:
        logging.error(f"Error in add_metadata: {str(e)}")
        logging.error(f"Error occurred on line: {e.__traceback__.tb_lineno}")
        raise

def main():
    try:
        logging.info("Processing data for cell_ID")
        
        pooled_data, intermediate_files = load_and_process_files(BASE_PATH, FILENAMES)
        
        if pooled_data.empty:
            logging.warning("No data to process for cell_ID. Exiting.")
            return

        logging.info(f"Pooled data shape: {pooled_data.shape}")
        
        pooled_data_with_metadata = add_metadata(BASE_PATH, pooled_data)

        output_file = os.path.join(BASE_PATH, 'cell_ID_pooled_median.txt')
        pooled_data_with_metadata.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Output saved to {output_file}")

        # Clean up intermediate files
        for f in intermediate_files:
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        logging.info("Intermediate files cleaned up")

    except Exception as e:
        logging.error(f"Error in main function: {str(e)}")
        logging.error(f"Error occurred on line: {e.__traceback__.tb_lineno}")
        raise

if __name__ == "__main__":
    main()