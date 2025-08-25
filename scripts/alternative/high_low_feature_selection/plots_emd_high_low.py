import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import gridspec

def load_and_process_data(file_path, concentrations=None):
    data = pd.read_csv(file_path, sep='\t')
    if concentrations:
        data = data[data['Concentration'].isin(concentrations)]
    return data

def create_emd_comparison_plot(treatment_data, control_data, feature_order, feature_positions, output_file, figsize=(28.8, 15), is_small=False):
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3], hspace=0.05)

    ax1 = plt.subplot(gs[0])  # Top subplot (log scale)
    ax2 = plt.subplot(gs[1])  # Bottom subplot (linear scale)

    dot_size = 0.2 if is_small else 6
    alpha = 0.6 if is_small else 0.5
    fontsize_factor = 0.8 if is_small else 1.5

    # Calculate control averages and IQR
    control_averages = control_data.groupby('Feature')['EMD_score'].mean()
    iqr_threshold = control_averages.quantile(0.75)
    inactivity_threshold = 0.05

    for feature in feature_order:
        # Plot treatment data
        treatment_feature_data = treatment_data[treatment_data['Feature'] == feature]
        x_treatment = [feature_positions[feature]] * len(treatment_feature_data)
        y_treatment = treatment_feature_data['EMD_score']
        
        mask_top_treatment = y_treatment > 7
        if np.any(mask_top_treatment):
            ax1.scatter(np.array(x_treatment)[mask_top_treatment], y_treatment[mask_top_treatment], c='gray', s=dot_size, alpha=alpha)
        
        mask_bottom_treatment = y_treatment <= 7
        ax2.scatter(np.array(x_treatment)[mask_bottom_treatment], y_treatment[mask_bottom_treatment], c='gray', s=dot_size, alpha=alpha)

        # Plot control data
        control_feature_data = control_data[control_data['Feature'] == feature]
        x_control = [feature_positions[feature]] * len(control_feature_data)
        y_control = control_feature_data['EMD_score']
        
        mask_top_control = y_control > 7
        if np.any(mask_top_control):
            ax1.scatter(np.array(x_control)[mask_top_control], y_control[mask_top_control], c='black', s=dot_size, alpha=alpha)
        
        mask_bottom_control = y_control <= 7
        ax2.scatter(np.array(x_control)[mask_bottom_control], y_control[mask_bottom_control], c='black', s=dot_size, alpha=alpha)

        # Plot control average
        avg_control = control_averages[feature]
        if avg_control > 7:
            ax1.scatter(feature_positions[feature], avg_control, c='black', s=dot_size*6, zorder=3)
        else:
            ax2.scatter(feature_positions[feature], avg_control, c='black', s=dot_size*6, zorder=3)

    # Plot IQR threshold
    if iqr_threshold > 7:
        ax1.axhline(y=iqr_threshold, color='#d7191c', linestyle='--', linewidth=1, zorder=2, dashes=(5, 5))
    else:
        ax2.axhline(y=iqr_threshold, color='#d7191c', linestyle='--', linewidth=1, zorder=2, dashes=(5, 5))

    # Plot inactivity threshold
    ax2.axhline(y=inactivity_threshold, color='#2c7bb6', linestyle='--', linewidth=1, zorder=2, dashes=(5, 5))

    ax1.set_yscale('log')
    ax1.set_ylim(7, max(treatment_data['EMD_score'].max(), control_data['EMD_score'].max()) * 1.1)
    ax1.set_xlim(-1, len(feature_order))
    ax1.set_xticks([])
    
    ax2.set_ylim(0, 7)
    ax2.set_xlim(-1, len(feature_order))

    # Remove only the top and right spines
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if is_small:
            ax.spines['left'].set_linewidth(0.5)
            ax.spines['bottom'].set_linewidth(0.5)

    if not is_small:
        ax2.set_xticks(range(len(feature_order)))
        ax2.set_xticklabels(feature_order, rotation=90, fontsize=4.2*fontsize_factor)
        ax2.tick_params(axis='x', which='major', pad=1)
        ax2.tick_params(axis='y', which='major', labelsize=14*fontsize_factor)
        ax1.tick_params(axis='y', which='major', labelsize=14*fontsize_factor)
        
        # Color the tick labels
        for tick, feature in zip(ax2.get_xticklabels(), feature_order):
            if control_averages[feature] > iqr_threshold:
                tick.set_color('#d7191c')
            elif all(treatment_data[treatment_data['Feature'] == feature]['EMD_score'] < inactivity_threshold):
                tick.set_color('#2c7bb6')
            else:
                tick.set_color('black')
    else:
        ax2.set_xticks([])
        ax2.set_xticklabels([])

    # Adjust subplot positions
    if is_small:
        plt.subplots_adjust(left=0.08, top=0.95, bottom=0.1)
    else:
        plt.subplots_adjust(left=0.06, top=0.95, bottom=0.15, right=0.99)

    # Add y-axis label
    if is_small:
        fig.text(0.02, 0.5, 'EMD score (linear/log scale)', fontsize=14*fontsize_factor, ha='left', va='center', rotation='vertical')
    else:
        fig.text(0.01, 0.5, 'EMD score (linear/log scale)', fontsize=14*fontsize_factor, ha='left', va='center', rotation='vertical')

    # Add 'Features' label
    if is_small:
        fig.text(0.5, 0.02, 'Features', fontsize=12*fontsize_factor, ha='center', va='center')
    else:
        ax2.set_xlabel('Features', fontsize=14*fontsize_factor, labelpad=10)

    fig.suptitle('EMD scores comparison', fontsize=16*fontsize_factor, y=0.98)

    # Add legend
    ax2.scatter([], [], c='gray', label='Treatment', s=dot_size*2)
    ax2.scatter([], [], c='black', label='Control', s=dot_size*2)
    ax2.scatter([], [], c='black', label='Control Average', s=dot_size*6)
    ax2.axhline(y=iqr_threshold, color='#d7191c', linestyle='--', linewidth=1, label='IQR Threshold', dashes=(5, 5))
    ax2.axhline(y=inactivity_threshold, color='#2c7bb6', linestyle='--', linewidth=1, label='Inactivity Threshold', dashes=(5, 5))
    ax2.legend(fontsize=10*fontsize_factor, loc='upper right')

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def select_features(treatment_data, control_data):
    control_averages = control_data.groupby('Feature')['EMD_score'].mean()
    iqr_threshold = control_averages.quantile(0.75)
    inactivity_threshold = 0.05

    selected_features = []
    for feature in control_averages.index:
        if control_averages[feature] <= iqr_threshold:
            feature_treatment_data = treatment_data[treatment_data['Feature'] == feature]
            if any(feature_treatment_data['EMD_score'] >= inactivity_threshold):
                selected_features.append(feature)

    return selected_features

# Define file paths and output directory
treatment_file = 'PATH/TO/feature_selection/results/EMD_scores/EMD_conc_2.5_97.5_well.txt'
control_file = 'PATH/TO//feature_selection/results/EMD_scores/EMD_c11_2.5_97.5_well.txt'
output_dir = 'PATH/TO/feature_selection_high_low/results/EMD_scores'
os.makedirs(output_dir, exist_ok=True)

# Load and process data
concentrations = ['2_vs_11', '3_vs_11', '4_vs_11', '5_vs_11', '6_vs_11', '7_vs_11', '8_vs_11', '9_vs_11', '10_vs_11']
treatment_data = load_and_process_data(treatment_file, concentrations)
control_data = load_and_process_data(control_file)

# Select features based on criteria
selected_features = select_features(treatment_data, control_data)

# Save selected features to a file
selected_features_file = os.path.join(output_dir, 'selected_features.txt')
with open(selected_features_file, 'w') as f:
    for feature in selected_features:
        f.write(f"{feature}\n")

# Combine data for feature ordering
combined_data = pd.concat([treatment_data, control_data])
feature_order = combined_data.groupby('Feature')['EMD_score'].median().sort_values().index
feature_positions = {feature: i for i, feature in enumerate(feature_order)}

# Create full-size plot
full_size_output = os.path.join(output_dir, 'EMD_comparison_scatter_plot_full_size.png')
create_emd_comparison_plot(treatment_data, control_data, feature_order, feature_positions, full_size_output)

# Create smaller plot for academic article
small_size_output = os.path.join(output_dir, 'EMD_comparison_scatter_plot_small_size.png')
create_emd_comparison_plot(treatment_data, control_data, feature_order, feature_positions, small_size_output, figsize=(8.64, 4.5), is_small=True)

print(f'Full-size EMD comparison scatter plot saved to: {full_size_output}')
print(f'Small-size EMD comparison scatter plot saved to: {small_size_output}')
print(f'List of selected features saved to: {selected_features_file}')
print(f'Number of selected features: {len(selected_features)}')