import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.pairwise import rbf_kernel
from scipy.spatial.distance import mahalanobis
from scipy.stats import wasserstein_distance
import matplotlib.pyplot as plt
import os
import glob
from sklearn.neighbors import KernelDensity
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['dim_reduction_visualization']['base_path']
subsets_min_count_filtered_dir_rel = config['dim_reduction_visualization']['subsets_min_count_filtered_dir']
mmd_emd_mah_output_dir_rel = config['dim_reduction_visualization']['mmd_emd_mah_output_dir']

# Construct full paths
input_folder = os.path.join(base_path, subsets_min_count_filtered_dir_rel)
output_folder = os.path.join(base_path, mmd_emd_mah_output_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

def replace_inf_and_large_values(df, max_value=1e12):
    df = df.copy()
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    for col in df.columns:
        df.loc[df[col] > max_value, col] = np.nan
        df.loc[df[col] < -max_value, col] = np.nan
    return df

def compute_mmd(X, Y, gamma=0.24):
    K_XX = rbf_kernel(X, X, gamma=gamma)
    K_YY = rbf_kernel(Y, Y, gamma=gamma)
    K_XY = rbf_kernel(X, Y, gamma=gamma)
    return np.mean(K_XX) + np.mean(K_YY) - 2 * np.mean(K_XY)

def compute_mahalanobis_distance(test_group, control_group):
    """
    Compute mean Mahalanobis distance between test group and control group (population 11)
    using control group as the reference distribution
    """
    mean_control = np.mean(control_group, axis=0)
    cov_control = np.cov(control_group, rowvar=False)
    
    try:
        inv_cov_control = np.linalg.inv(cov_control)
    except np.linalg.LinAlgError:
        cov_control += np.eye(cov_control.shape[0]) * 1e-6
        inv_cov_control = np.linalg.inv(cov_control)
    
    distances = []
    for x in test_group:
        try:
            dist = mahalanobis(x, mean_control, inv_cov_control)
            distances.append(dist)
        except ValueError:
            continue
    
    return np.mean(distances) if distances else np.nan

def compute_wasserstein(X, Y):
    """
    Compute mean Wasserstein distance (EMD) across all features
    """
    distances = []
    for feature_idx in range(X.shape[1]):
        dist = wasserstein_distance(X[:, feature_idx], Y[:, feature_idx])
        distances.append(dist)
    return np.mean(distances)

def process_file(file_path):
    data = pd.read_csv(file_path, delimiter='\t')
    
    exclude_cols = [
        'Concentration', 'counts_Cells', 'counts_Cytoplasm', 'counts_FilteredNuclei', 
        'Metadata_Well', 'Metadata_Day', 'Metadata_Biorep', 'Tech_replica', 'Day_Well_BR', 'cell_ID'
    ]
    
    feature_cols = [col for col in data.columns if col not in exclude_cols]
    
    data[feature_cols] = replace_inf_and_large_values(data[feature_cols])
    
    for concentration in data['Concentration'].unique():
        group_mask = data['Concentration'] == concentration
        for col in feature_cols:
            mean_value = data.loc[group_mask, col].mean()
            data.loc[group_mask, col] = data.loc[group_mask, col].fillna(mean_value)
    
    data[feature_cols] = replace_inf_and_large_values(data[feature_cols])
    
    nan_features = data[feature_cols].columns[data[feature_cols].isna().any()].tolist()
    if nan_features:
        print(f"Features with NaN values: {nan_features}")
        feature_cols = [col for col in feature_cols if col not in nan_features]
        if not feature_cols:
            raise ValueError("All features contain NaN values after processing.")
    
    data = data[['Concentration'] + feature_cols]
    
    scaler = MinMaxScaler()
    data[feature_cols] = scaler.fit_transform(data[feature_cols])
    
    comparisons = [(11, 10), (11, 9), (11, 8), (11, 7), (11, 6), (11, 5), (11, 4), (11, 3), (11, 2)]
    results = []
    
    for (a, b) in comparisons:
        group_a = data[data['Concentration'] == a][feature_cols].values
        group_b = data[data['Concentration'] == b][feature_cols].values
        
        results.append({
            'Group_A': a,
            'Group_B': b,
            'MMD_Score': compute_mmd(group_a, group_b),
            'Mahalanobis_Score': compute_mahalanobis_distance(group_b, group_a),
            'EMD_Score': compute_wasserstein(group_a, group_b)
        })
    
    return pd.DataFrame(results)

def create_metric_plot(data, metric, title, output_path, color='#80c0c0'):
    plt.figure(figsize=(12, 8))
    plt.plot(range(len(data)), data[metric], marker='o', color=color, linewidth=2)
    plt.xticks(range(len(data)), [f'{a} vs {b}' for a, b in zip(data['Group_A'], data['Group_B'])])
    plt.xlabel('Comparisons')
    plt.ylabel(metric.replace('_', ' '))
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_combined_plot(data, metrics, title, output_path):
    plt.figure(figsize=(15, 10))
    
    # Create a plot with multiple y-axes
    fig, axes = plt.subplots(figsize=(15, 10))
    
    colors = ['#80c0c0', '#333333', '#777777']
    
    for idx, (metric, color) in enumerate(zip(metrics, colors)):
        if idx == 0:
            ax = axes
        else:
            ax = axes.twinx()
            ax.spines['right'].set_position(('outward', 60 * (idx - 1)))
        
        line = ax.plot(range(len(data)), data[metric], 
                      marker='o', label=metric.replace('_', ' '), 
                      color=color, linewidth=2)
        ax.set_ylabel(metric.replace('_', ' '), color=color)
        ax.tick_params(axis='y', labelcolor=color)
    
    plt.xlabel('Comparisons')
    plt.xticks(range(len(data)), [f'{a} vs {b}' for a, b in zip(data['Group_A'], data['Group_B'])])
    plt.title(title)
    
    # Combine legends
    lines = []
    labels = []
    for ax in fig.axes:
        axline, axlabel = ax.get_legend_handles_labels()
        lines.extend(axline)
        labels.extend(axlabel)
    
    plt.legend(lines, labels, loc='center left', bbox_to_anchor=(1.2, 0.5))
    
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

# Main execution
colors = ['#80c0c0', '#333333', '#777777']

# Full color palette for all_conditions plots - keeping 4 colors as requested
all_conditions_colors = ['#80c0c0', '#008080', '#777777', '#333333']

day_conversion = {
    'D1': 'Day 1',
    'D5': 'Day 5',
    'D7': 'Day 7',
    'D9': 'Day 9'
}

metrics = ['MMD_Score', 'Mahalanobis_Score', 'EMD_Score']

all_results = []
for file_path in glob.glob(os.path.join(input_folder, '*.txt')):
    file_name = os.path.basename(file_path)
    condition = file_name.split('_')[2]
    
    try:
        print(f"Processing file: {file_name}")
        results_df = process_file(file_path)
        results_df['Condition'] = condition
        all_results.append(results_df)
        
        # Plot individual metrics
        for metric, color in zip(metrics, colors):
            title = f'{metric.replace("_", " ")} - {day_conversion.get(condition, condition)}'
            output_path = os.path.join(output_folder, f'{metric.lower()}_plot_{condition}.png')
            create_metric_plot(results_df, metric, title, output_path, color)
        
        # Plot combined metrics
        title = f'All Metrics - {day_conversion.get(condition, condition)}'
        output_path = os.path.join(output_folder, f'combined_plot_{condition}.png')
        create_combined_plot(results_df, metrics, title, output_path)
        
    except Exception as e:
        print(f"Error processing file {file_name}: {str(e)}")

if not all_results:
    print("No files were successfully processed. Check the input data and error messages above.")
else:
    consolidated_results = pd.concat(all_results, ignore_index=True)
    consolidated_results.to_csv(os.path.join(output_folder, 'consolidated_results.txt'), sep='\t', index=False)

    # Plot all conditions for each metric
    for metric, color in zip(metrics, colors):
        plt.figure(figsize=(12, 8))
        conditions = sorted(consolidated_results['Condition'].unique())

        for i, condition in enumerate(conditions):
            condition_data = consolidated_results[consolidated_results['Condition'] == condition]
            plt.plot(range(len(condition_data)), condition_data[metric], 
                    marker='o', label=day_conversion.get(condition, condition), 
                    color=all_conditions_colors[i % len(all_conditions_colors)], linewidth=4)

        plt.xticks(range(len(condition_data)), 
                  [f'{a} vs {b}' for a, b in zip(condition_data['Group_A'], condition_data['Group_B'])])
        plt.xlabel('Comparisons')
        plt.ylabel(metric.replace('_', ' '))
        plt.title(f'{metric.replace("_", " ")} - All Conditions')
        
        legend = plt.legend(title="Days:", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.setp(legend.get_title(), fontweight='bold')
        
        plt.tight_layout()
        plt.grid(True)
        plt.savefig(os.path.join(output_folder, f'{metric.lower()}_all_conditions.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()

print(f"Processing complete. Results saved in {output_folder}")