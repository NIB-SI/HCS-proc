import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import os
from multiprocessing import Pool

def replace_inf_and_large_values(df, max_value=1e12):
    """Replace infinite and very large values with NaN"""
    df = df.copy()
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    for col in df.columns:
        if df[col].dtype in ['float64', 'float32', 'int64', 'int32']:
            df.loc[df[col] > max_value, col] = np.nan
            df.loc[df[col] < -max_value, col] = np.nan
    return df

def significant_diff_from_others(values, threshold=0.3):
    """Check if any mean value is significantly different from others."""
    for value in values:
        for other in values:
            if other != value and other != 0:
                if abs(value - other) / other > threshold:
                    return True
    return False

def plot_column(col, df_path, plots_dir):
    """Plot violin plots for a single column across all concentrations for each day"""
    df = pd.read_csv(df_path, sep='\t')
    
    # Handle infinite and large values
    df = replace_inf_and_large_values(df)
    
    # Drop rows where the current column has NaN values
    df.dropna(subset=[col], inplace=True)
    
    if df.empty:
        print(f"No valid data for column {col}, skipping...")
        return
    
    # Create directory for this column
    col_dir = os.path.join(plots_dir, col)
    os.makedirs(col_dir, exist_ok=True)
    
    # Get unique days from Metadata_Day
    available_days = sorted(df['Metadata_Day'].unique())
    
    # Plot for each day
    for day in available_days:
        day_data = df[df['Metadata_Day'] == day]
        
        if day_data.empty:
            continue
            
        plt.figure(figsize=(12, 6))
        
        # Get available concentrations for this day, sorted
        available_concentrations = sorted(day_data['Concentration'].unique())
        
        # Create violin plot
        ax = sns.violinplot(x='Concentration', y=col, data=day_data, 
                           inner=None, order=available_concentrations, color='#82cafc')
        
        # Calculate means and medians for each concentration
        group_stats = day_data.groupby('Concentration')[col].agg(['median', 'mean'])
        
        # Plot means and medians
        for i, concentration in enumerate(available_concentrations):
            if concentration in group_stats.index:
                # Median (yellow X)
                median_val = group_stats.loc[concentration]['median']
                ax.scatter(i, median_val, color='#fddc5c', edgecolor='black', 
                          marker='X', s=50, zorder=3)
                # Mean (red circle)
                mean_val = group_stats.loc[concentration]['mean']
                ax.scatter(i, mean_val, color='#fc5a50', edgecolor='black', zorder=3)
        
        # Annotate 'n' values below the violins
        for i, concentration in enumerate(available_concentrations):
            n = day_data[day_data['Concentration'] == concentration][col].dropna().shape[0]
            ax.text(i, ax.get_ylim()[0], f'n={n}', horizontalalignment='center', 
                   size='x-small', color='black', weight='normal')
        
        # Reset position of 'Concentration' label
        ax.xaxis.labelpad = 10
        
        # Perform Kruskal-Wallis test and annotate significant differences
        groups = [group[col].dropna() for name, group in day_data.groupby('Concentration')]
        title_text = f'{col} - Day {day}'
        
        if len(groups) > 1 and any(len(group) > 0 and group.nunique() > 1 for group in groups):
            try:
                # Only perform test if we have valid groups
                valid_groups = [group for group in groups if len(group) > 0 and group.nunique() > 1]
                if len(valid_groups) > 1:
                    h, p = kruskal(*valid_groups)
                    if p < 0.05:
                        title_text += f' (p<0.05)'
                    
                    # Check for significant differences in means
                    if significant_diff_from_others(group_stats['mean'].values, threshold=0.3):
                        title_text += ' (M)'
                        # Mark concentrations with significant mean differences
                        for i, (concentration, mean) in enumerate(zip(available_concentrations, group_stats['mean'].values)):
                            for other_mean in group_stats['mean'].values:
                                if other_mean != mean and other_mean != 0:
                                    if abs(mean - other_mean) / abs(other_mean) > 0.3:
                                        ax.text(i, ax.get_ylim()[1] * 0.95, '*', 
                                               horizontalalignment='center', color='black', weight='bold')
                                        break
            except (ValueError, ZeroDivisionError) as e:
                print(f"Statistical test failed for {col}, Day {day}: {e}")
        
        plt.title(title_text, fontsize=14, y=1.05)
        plt.xlabel('Concentration')
        plt.ylabel(col)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(col_dir, f'cell_new_{col}_Day_{day}.png'), dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    file_path = '/PATH/TO/feature_selection/results/trimmed_2/clean_trimmed_features_all_days_trimmed_trimmed_features_cid.txt'
    plots_dir = '/PATH/TO/violin_plots/results'
    os.makedirs(plots_dir, exist_ok=True)
    
    # Read the file to get column names
    df = pd.read_csv(file_path, sep='\t')
    
    # Define columns to exclude (metadata and count columns)
    exclude_cols = [
        'Concentration', 'counts_Cells', 'counts_Cytoplasm', 'counts_FilteredNuclei', 
        'Metadata_Well', 'Metadata_Day', 'Metadata_Biorep', 'Tech_replica', 'Day_Well_BR', 'cell_ID'
    ]
    
    # Get feature columns (everything except excluded columns)
    feature_columns = [col for col in df.columns if col not in exclude_cols]
    
    print(f"Found {len(feature_columns)} feature columns to plot")
    print(f"Available days: {sorted(df['Metadata_Day'].unique())}")
    print(f"Available concentrations: {sorted(df['Concentration'].unique())}")
    
    # Use multiprocessing to plot columns in parallel
    num_processes = 20
    with Pool(num_processes) as pool:
        pool.starmap(plot_column, [(col, file_path, plots_dir) for col in feature_columns])
    
    print(f"Plotting complete. Results saved in {plots_dir}")