import os
import pandas as pd
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['categorization']['base_path']
final_features_file = os.path.join(base_path, config['categorization']['final_features_file'])
output_dir = os.path.join(base_path, config['categorization']['output_dir'])
categories_file = config['categorization']['feature_categories_file']

# Table layout - the display order of the rows and columns, and how the categories and
# organelles recorded in the lookup map onto them
category_columns = [c.strip() for c in config['categorization']['category_columns'].split(',')]
organelle_rows = [r.strip() for r in config['categorization']['organelle_rows'].split(',')]
category_labels = dict(item.split(':') for item in config['categorization']['category_labels'].split(','))
organelle_labels = dict(item.split(':') for item in config['categorization']['organelle_labels'].split(','))

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Step 1: Read the selected features and the universal categorization
with open(final_features_file) as fh:
    selected = [line.strip() for line in fh if line.strip()]

categories = pd.read_csv(categories_file, sep='\t', comment='#', dtype=str).fillna('')
lookup = {r.feature: (r.category, r.organelle) for r in categories.itertuples()}
print(f"{len(selected)} selected features, {len(lookup)} in {categories_file}")

# Step 2: Look every selected feature up, reporting anything the table cannot resolve
counts = pd.DataFrame(0, index=organelle_rows, columns=category_columns)
unknown, undecided = [], []

for feature in selected:
    if feature not in lookup:
        unknown.append(feature)
        continue
    category, organelle = lookup[feature]
    column = category_labels.get(category)
    row = organelle_labels.get(organelle)
    if not column or not row:
        undecided.append(feature)
        continue
    counts.at[row, column] += 1

if unknown:
    print(f"\n{len(unknown)} feature(s) missing from {categories_file} - add them with "
          f"build_feature_categories.py, then fill in the assignment:")
    for feature in unknown:
        print(f"  {feature}")
if undecided:
    print(f"\n{len(undecided)} feature(s) have no category or organelle assigned yet:")
    for feature in undecided:
        print(f"  {feature}")

counted = int(counts.values.sum())
if counted != len(selected):
    print(f"\nWARNING: {counted} of {len(selected)} selected features were counted.")

# Step 3: Save the cross-tab. Zeros are written as blanks so the table reads the way it
# is plotted - an absent combination is not a measured zero
output_file = os.path.join(output_dir, 'types.txt')
counts.replace(0, '').to_csv(output_file, sep='\t')

print(f"\nFeature counts ({counted} features):")
print(counts.to_string())
print("\nCounts per category:")
print(counts.sum(axis=0).to_string())
print(f"\nCategorization table saved to {output_file}")
