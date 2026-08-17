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

# Step 2: Look every selected feature up. Anything that cannot be placed is collected
# with the reason, so a gap in the lookup is never mistaken for a genuine zero.
counts = pd.DataFrame(0, index=organelle_rows, columns=category_columns)
unplaced = []

for feature in selected:
    if feature not in lookup:
        unplaced.append((feature, f"not in {os.path.basename(categories_file)}"))
        continue
    category, organelle = lookup[feature]
    column = category_labels.get(category)
    row = organelle_labels.get(organelle)
    reasons = []
    if not category:
        reasons.append("category blank")
    elif not column:
        reasons.append(f"category '{category}' is not in category_labels")
    if not organelle:
        reasons.append("organelle blank")
    elif not row:
        reasons.append(f"organelle '{organelle}' is not in organelle_labels")
    if reasons:
        unplaced.append((feature, "; ".join(reasons)))
        continue
    counts.at[row, column] += 1

# Step 3: Save the cross-tab. Zeros are written as blanks so the table reads the way it
# is plotted - an absent combination is not a measured zero
counted = int(counts.values.sum())
output_file = os.path.join(output_dir, 'types.txt')
counts.astype(object).where(counts != 0, '').to_csv(output_file, sep='\t')

print(f"\nFeature counts ({counted} features):")
print(counts.to_string())
print("\nCounts per category:")
print(counts.sum(axis=0).to_string())
print(f"\nCategorization table saved to {output_file}")

# Step 4: Report the gaps last, so they are the final thing on screen
if unplaced:
    print(f"\n{'=' * 70}")
    print(f"ACTION NEEDED: {len(unplaced)} of {len(selected)} selected features were not "
          f"counted.\nFix them in {categories_file} - run build_feature_categories.py to "
          f"append\nany that are missing, then fill in the blanks by hand.")
    print('=' * 70)
    for feature, reason in unplaced:
        print(f"  {feature}\n      {reason}")
else:
    print(f"\nAll {len(selected)} selected features were categorized.")
