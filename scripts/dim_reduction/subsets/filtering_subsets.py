import pandas as pd
import numpy as np
import os
import shutil
from pathlib import Path
import time
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_subsets']['base_path']
subsets_all_days_fs_dir_rel = config['dim_reduction_subsets']['subsets_all_days_fs_dir']
subsets_filtered_dir_rel = config['dim_reduction_subsets']['subsets_filtered_dir']

# Construct full paths
input_dir = os.path.join(base_path, subsets_all_days_fs_dir_rel)
output_dir = os.path.join(base_path, subsets_filtered_dir_rel)

# Define the metadata columns that should not be filtered
metadata_columns = [
    'Concentration', 'counts_Cells', 'counts_Cytoplasm', 'counts_FilteredNuclei',
    'Metadata_Well', 'Metadata_Day', 'Metadata_Biorep', 'Tech_replica', 'Day_Well_BR', 'cell_ID'
]

# Define specific features to exclude regardless of other criteria
features_to_exclude = [
    'rp_norm_AreaShape_Compactness_RelateLysoCell',
    'rp_norm_Texture_AngularSecondMoment_GrayLys_3_00_256_RelateLysoCell'
]

# Create the report file
os.makedirs(output_dir, exist_ok=True)
report_path = os.path.join(output_dir, 'feature_filtering_report.txt')

# Initialize counters for summary statistics as global variables
total_files_processed = 0
total_features_removed = 0
total_nan_removed = 0
total_zero_var_removed = 0
total_obj_number_removed = 0
total_excluded_removed = 0

with open(report_path, 'w') as report_file:
    report_file.write("Feature Filtering Report\n")
    report_file.write("======================\n")
    report_file.write(f"Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    report_file.write(f"Input Directory: {input_dir}\n")
    report_file.write(f"Output Directory: {output_dir}\n\n")
    report_file.write("Filtering Criteria:\n")
    report_file.write("1. Features with ≥30% NaN values\n")
    report_file.write("2. Features with zero variance\n")
    report_file.write("3. Features containing 'Number_Object_Number' in their name\n")
    report_file.write("4. Explicitly excluded features\n\n")
    
    # Write the list of explicitly excluded features to the report
    report_file.write("Explicitly excluded features:\n")
    for feature in features_to_exclude:
        report_file.write(f"  - {feature}\n")
    report_file.write("\n")
    
    # Function to process a single file
    def process_file(file_path, rel_path):
        global total_files_processed, total_features_removed, total_nan_removed
        global total_zero_var_removed, total_obj_number_removed, total_excluded_removed
        
        # Read the file
        data = pd.read_csv(file_path, sep='\t')
        original_feature_count = data.shape[1] - len(metadata_columns)
        
        # Separate features and metadata
        features = data.drop(columns=metadata_columns, errors='ignore')
        metadata = data[metadata_columns].copy()
        
        # Initialize lists to track removed features
        nan_features = []
        zero_var_features = []
        obj_number_features = []
        excluded_features = [f for f in features_to_exclude if f in features.columns]
        
        # Check for features with ≥30% NaN values
        nan_counts = features.isna().sum()
        nan_percentages = (nan_counts / len(features)) * 100
        nan_features = list(nan_percentages[nan_percentages >= 30.0].index)
        
        # Replace inf values with NaN for variance calculation
        features_inf_replaced = features.replace([np.inf, -np.inf], np.nan)
        
        # Check for features with zero variance
        variances = features_inf_replaced.var(skipna=True)
        zero_var_features = list(variances[variances == 0].index)
        
        # Check for features containing "Number_Object_Number"
        obj_number_features = [col for col in features.columns if "Number_Object_Number" in col]
        
        # Combine all features to remove
        features_to_remove = list(set(nan_features + zero_var_features + obj_number_features + excluded_features))
        
        # Update counters
        total_features_removed += len(features_to_remove)
        total_nan_removed += len(nan_features)
        total_zero_var_removed += len(zero_var_features)
        total_obj_number_removed += len(obj_number_features)
        total_excluded_removed += len(excluded_features)
        
        # Remove the features
        filtered_data = data.drop(columns=features_to_remove, errors='ignore')
        
        # Create output directory if it doesn't exist
        output_file_dir = os.path.dirname(os.path.join(output_dir, rel_path))
        os.makedirs(output_file_dir, exist_ok=True)
        
        # Save the filtered data
        output_file_path = os.path.join(output_dir, rel_path)
        filtered_data.to_csv(output_file_path, sep='\t', index=False)
        
        # Write to the report
        report_file.write(f"File: {rel_path}\n")
        report_file.write(f"  Original feature count: {original_feature_count}\n")
        report_file.write(f"  Removed feature count: {len(features_to_remove)}\n")
        
        if nan_features:
            report_file.write(f"  Features removed due to ≥30% NaN values ({len(nan_features)}):\n")
            for feature in sorted(nan_features):
                report_file.write(f"    - {feature}: {nan_percentages[feature]:.2f}%\n")
        
        if zero_var_features:
            report_file.write(f"  Features removed due to zero variance ({len(zero_var_features)}):\n")
            for feature in sorted(zero_var_features):
                report_file.write(f"    - {feature}\n")
        
        if obj_number_features:
            report_file.write(f"  Features removed containing 'Number_Object_Number' ({len(obj_number_features)}):\n")
            for feature in sorted(obj_number_features):
                report_file.write(f"    - {feature}\n")
        
        if excluded_features:
            report_file.write(f"  Explicitly excluded features ({len(excluded_features)}):\n")
            for feature in sorted(excluded_features):
                report_file.write(f"    - {feature}\n")
        
        report_file.write(f"  Remaining feature count: {filtered_data.shape[1] - len(metadata_columns)}\n\n")
        
        return len(features_to_remove)
    
    # Walk through the directory structure
    for root, dirs, files in os.walk(input_dir):
        # Get the relative path to maintain directory structure
        rel_root = os.path.relpath(root, input_dir)
        
        # Process each .txt file
        for file in files:
            if file.endswith('.txt'):
                file_path = os.path.join(root, file)
                rel_path = os.path.join(rel_root, file) if rel_root != '.' else file
                
                try:
                    features_removed = process_file(file_path, rel_path)
                    total_files_processed += 1
                    print(f"Processed: {rel_path} - Removed {features_removed} features")
                except Exception as e:
                    print(f"Error processing {rel_path}: {str(e)}")
                    report_file.write(f"Error processing {rel_path}: {str(e)}\n\n")
    
    # Write summary statistics
    report_file.write("\nSummary Statistics\n")
    report_file.write("=================\n")
    report_file.write(f"Total files processed: {total_files_processed}\n")
    report_file.write(f"Total features removed: {total_features_removed}\n")
    report_file.write(f"  - Due to ≥30% NaN values: {total_nan_removed}\n")
    report_file.write(f"  - Due to zero variance: {total_zero_var_removed}\n")
    report_file.write(f"  - Containing 'Number_Object_Number': {total_obj_number_removed}\n")
    report_file.write(f"  - Explicitly excluded: {total_excluded_removed}\n")

print(f"Feature filtering complete! Report saved to {report_path}")