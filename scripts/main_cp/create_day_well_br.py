import os
import logging
import configparser
import pandas as pd
import shutil

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

def add_day_well_br(df: pd.DataFrame) -> pd.DataFrame:
    if all(col in df.columns for col in ['Metadata_Day', 'Metadata_Well', 'Metadata_Biorep']):
        df['Day_Well_BR'] = (df['Metadata_Day'].astype(str) + '_' + 
                             df['Metadata_Well'].astype(str) + '_' + 
                             df['Metadata_Biorep'].astype(str))
    else:
        logging.warning(f"Required columns for Day_Well_BR not found in dataframe")
    return df

def process_file(filepath: str):
    logging.info(f"Processing file: {filepath}")
    
    # Read the file
    df = pd.read_csv(filepath, sep='\t')
    
    # Add Day_Well_BR column
    df = add_day_well_br(df)
    
    # Create a temporary file with the updated data
    temp_filepath = filepath + '.temp'
    df.to_csv(temp_filepath, sep='\t', index=False)
    
    logging.info(f"Updated file saved temporarily: {temp_filepath}")
    
    return temp_filepath

def main():
    temp_files = []
    try:
        for filename in FILENAMES:
            filepath = os.path.join(BASE_PATH, filename)
            
            if not os.path.exists(filepath):
                logging.warning(f"File not found: {filepath}")
                continue
            
            temp_filepath = process_file(filepath)
            temp_files.append((filepath, temp_filepath))

        # Replace original files with updated ones
        for original_file, temp_file in temp_files:
            shutil.move(temp_file, original_file)
            logging.info(f"Replaced {original_file} with updated version")

        logging.info("All files processed and replaced successfully")

    except Exception as e:
        logging.error(f"Error in main function: {str(e)}")
        logging.error(f"Error occurred on line: {e.__traceback__.tb_lineno}")
        # Clean up any temporary files if an error occurred
        for _, temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        raise

if __name__ == "__main__":
    main()